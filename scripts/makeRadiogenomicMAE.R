library(dplyr)
library(MuData)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(readxl)
library(rhdf5)
library(stringr)
library(tools)
library(yaml)

#' Function to make the radiomic Summarized Experiment object containing assays
#' for each image filter, image metadata as the coldata, and the pyradiomic
#' configuration as metadata.
#'
#' @param allRadiomicFeatures A dataframe with samples as rows, image metadata as first columns, remaining columns pyradiomic features.
#' @param pyradiomicsConfig A list read in from the yaml file used for radiomic feature extraction.
#' @param findFeature A string, a feature that has been included in the extraction, used to find all the filter types.
#' @examples
#' makeRadiomicSEO(nsclcRadiomicFeatures, findFeature = 'firstorder_10Percentile, pyradiomicsConfig = 'pyradSettings')

makeRadiomicSEO <- function(allRadiomicFeatures,
                            pyradiomicsConfig = NULL,
                            findFeature = "firstorder_10Percentile") { # nolint

    # Add check that the columns are the radiomic features and rows are the patient IDs

    # Get list of filter types included in this batch of radiomics
    filterList <- str_subset(names(allRadiomicFeatures), findFeature)
    # Get just the names of each of the filters
    filterList <- str_remove_all(filterList, findFeature)

    # Extract the shape features from the total list - these only get measured once
    shapeFeatures <- select(allRadiomicFeatures, contains("shape", ignore.case = TRUE))
    # Confirm if there are shape features included in the radiomic set
    if (dim(shapeFeatures)[2] == 0) {
        hasShapeFlag <- FALSE
    } else {
        hasShapeFlag <- TRUE
        shapeFeatureNames <- names(shapeFeatures)
        # names(shapeFeatures) <- str_remove_all(shapeFeatureNames, "original_")
    }

    # Drop shape features from total feature list
    noShapeRadiomicFeatures <- allRadiomicFeatures[, !names(allRadiomicFeatures) %in% shapeFeatureNames]

    # Extract image metadata from front of dataframe
    diagnosticInfoIndexStart <- grep("diagnostics_", colnames(noShapeRadiomicFeatures))[1]
    patientIDCol <- grep("patient_ID", colnames(noShapeRadiomicFeatures))
    imageMetadata <- noShapeRadiomicFeatures[patientIDCol:diagnosticInfoIndexStart - 1]
    # Combine the patient ID and ROI name to make unique rownames
    rownames(imageMetadata) <- paste(imageMetadata$patient_ID, imageMetadata$roi, sep = "_")

    # Initialize the list of filter assays for storage in loop
    filterAssays <- list()
    # Create an assay for each of the filters to be joined together in an experiment
    for (filter in filterList) {
        # Get features for specific filter
        filterRadiomicFeatures <- select(noShapeRadiomicFeatures,
                                         starts_with(filter, ignore.case = TRUE))
        # Have to remove filter prefix from feature names so all assay row/col names match
        names(filterRadiomicFeatures) <- str_remove_all(names(filterRadiomicFeatures), filter)
        # Append shape features to front of the radiomic feature dataframe for this filter
        if (hasShapeFlag == TRUE){
            filterRadiomicFeatures <- cbind(shapeFeatures, filterRadiomicFeatures)
        }

        # convert feature dataframe to a matrix and transpose so columns are patients and rows are features
        matFilterRadiomicFeatures <- t(data.matrix(filterRadiomicFeatures))
        # set column names to patient IDs
        colnames(matFilterRadiomicFeatures) <- rownames(imageMetadata)

        # Add this set of filter features to the list
        filterAssays <- append(filterAssays, list(matFilterRadiomicFeatures))
    }
    # Assign names of each filter to the list
    names(filterAssays) <- filterList

    # Construct summarized experiment
    radiomicFeaturesSEO <- SummarizedExperiment(assays = filterAssays,
                                                colData = imageMetadata,
                                                metadata = pyradiomicsConfig)

    return(radiomicFeaturesSEO)
}



#' Function to load in the clinical, genomic, or radiomic data for the radiogenomic MAE object
#'
#' @param dataFilePath A string path to the file to load.
loadRadiogenomicDataFile <- function(dataFilePath) {
    dataFileType = file_ext(dataFilePath)
    if (length(dataFileType) == 0) {
        # No file passed, return 0
        return(NULL)
    }
    if (dataFileType == "csv") {
        loadedDataframe <- read.csv(dataFilePath, header = TRUE, sep = ",", check.names = FALSE)
    } else if (dataFileType == "xlsx") {
        loadedDataframe <- read_excel(dataFilePath)
    } else {
        stop("Radiogenomic data file must be a .csv or .xlsx.")
    }

    return(loadedDataframe)
}


#' Function to construct a radiogenomic MultiAssayExperiment object.
#' Also built to build just a radiomic one with no genomic data
#'
#' @param clinicalDataFilePath A string of the path to the clinical data file (csv or xlsx)
#' @param radiomicDataFilePath A string of the path to the radiomic features file (csv or xlsx)
#' @param pyradiomicsConfigFile A string of the path to the configuration file for the PyRadiomics radiomic feature extraction (yaml)
#' @param findFeature A string, a feature that has been included in the extraction, used to find all the filter types.
#' @param genomicDataFilePath A string of the path to the genomic feature file (csv or xlsx). If not provided, radiomic MAE constructed.
#' @param outputFileName A string, file name to save the MAE to
#' @param clinicalPatIDCol A string, column name of the patient IDs in the clinical data file
#' @param radiomicPatIDCol A string, column name of the patient IDs in the radiomic feature file
makeRadiogenomicMAE <- function(clinicalDataFilePath,
                                radiomicDataFilePath,
                                pyradiomicsConfigFile,
                                findFeature = "firstorder_10Percentile",
                                genomicDataFilePath = NULL,
                                outputFileName = "outputMAE.rds",
                                clinicalPatIDCol = "Case ID",
                                radiomicPatIDCol = "patient_ID") {
    # Load everything
    clinicalDataframe = loadRadiogenomicDataFile(clinicalDataFilePath)
    radiomicDataframe = loadRadiogenomicDataFile(radiomicDataFilePath)
    genomicDataframe = loadRadiogenomicDataFile(genomicDataFilePath)

    if (is.null(genomicDataframe == 0)) {
        noGene = TRUE
        print("No genomic data passed. MAE will only contain radiomic data.")
    } else {
        noGene = FALSE
    }

    if (file_ext(pyradiomicsConfigFile) != "yaml") {
        stop("Radiomic configuration file must be of type yaml")
    } else {
        pyradiomicsConfig <- read_yaml(pyradiomicsConfigFile)
    }

    # Cross reference patient IDs
    clinicalPatIDs <- clinicalDataframe[[clinicalPatIDCol]]
    radiomicPatIDs <- radiomicDataframe[[radiomicPatIDCol]]
    # Get intersection of patients in clinical and radiomic data
    selectPatIDs <- unique(intersect(clinicalPatIDs, radiomicPatIDs))

    if (noGene != FALSE) {
        # ASSUMES THAT PATIENT ID IS A COLUMN IN RADIOMICS AND PATIENT IDS ARE COLUMN NAMES IN GENOMICS
        genePatientMatchesIdx <- match(selectPatIDs, names(genomicDataframe))
        genePatientMatchesIdx <- genePatientMatchesIdx[!is.na(genePatientMatchesIdx)]
        # Patient IDs in clinical, genomic and radiomic data
        selectPatIDs <- names(genomicDataframe)[genePatientMatchesIdx]

        selectPatGenomicData <- genomicDataframe[, c(names(genomicDataframe)[1], selectPatIDs)]

        # make this into a function at some point
        matSelectPatGenomicData <- data.matrix(selectPatGenomicData[, -(1)])
        rownames(matSelectPatGenomicData) <- selectPatGenomicData[, 1][[1]]
        genomicSEO <- SummarizedExperiment(assays = SimpleList(matSelectPatGenomicData),
                                           colData = selectPatIDs)
    }

    selectPatClinicalData <- filter(clinicalDataframe, clinicalDataframe[[clinicalPatIDCol]] %in% selectPatIDs)
    selectPatClinicalData <- as.data.frame(selectPatClinicalData)
    rownames(selectPatClinicalData) <- selectPatClinicalData[[clinicalPatIDCol]]

    selectPatRadiomicData <- filter(radiomicDataframe, radiomicDataframe[[radiomicPatIDCol]] %in% selectPatIDs)
    radiomicSEO <- makeRadiomicSEO(selectPatRadiomicData, pyradiomicsConfig = pyradiomicsConfig)

    # Construct sampleMap for MAE
    radiomicColname <- colnames(radiomicSEO)
    radiomicPrimary <- selectPatRadiomicData[[radiomicPatIDCol]]
    radiomicAssay <- rep(factor("radiomics"), each = length(radiomicPrimary))
    sampleMap <- data.frame(assay = radiomicAssay,
                            primary = radiomicPrimary,
                            colname = radiomicColname)

    if (noGene != FALSE) {
        experimentList <- ExperimentList(list(radiomics = radiomicSEO,
                                              rnaseq = genomicSEO))

        genomicColname <- colnames(genomicSEO)
        genomicPrimary <- selectPatIDs
        genomicAssay <- rep(factor("rnaseq"), each = length(genomicPrimary))
        genomicSampleMap <- data.frame(assay = genomicAssay,
                                       primary = genomicPrimary,
                                       colname = genomicColname)
        sampleMap <- rbind(sampleMap, genomicSampleMap)

    } else {
        experimentList <- ExperimentList(list(radiomics = radiomicSEO))
    }

    # preppedMAE <- prepMultiAssay(experimentList, selectPatClinicalData, sampleMap)

    # return(preppedMAE)

    radiogenomicMAE <- MultiAssayExperiment(experiments = experimentList,
                                            colData = selectPatClinicalData,
                                            sampleMap = sampleMap)

    outputExt <- file_ext(outputFileName)
    if (outputExt == "rds") {
        saveRDS(radiogenomicMAE, outputFileName)
    } else if (outputExt == "h5mu") {
        writeH5MU(radiogenomicMAE, outputFileName)
    }

    return(radiogenomicMAE)
}

# -- Read in Snakemake parameters
clinicalDataFilePath <- snakemake@input$clinical
radiomicDataFilePath <- snakemake@input$radiomic
pyradiomicsConfigFile <- snakemake@input$pyrad

findFeature <- snakemake@params$findFeature
clinicalPatIDCol <- snakemake@params$clinicalPatIDCol
radiomicPatIDCol <- snakemake@params$radiomicPatIDCol
outputFileName <- snakemake@params$outputFileName

# Call function
makeRadiogenomicMAE(clinicalDataFilePath,
                    radiomicDataFilePath,
                    pyradiomicsConfigFile,
                    findFeature,
                    outputFileName = outputFileName,
                    clinicalPatIDCol = clinicalPatIDCol,
                    radiomicPatIDCol = radiomicPatIDCol)
