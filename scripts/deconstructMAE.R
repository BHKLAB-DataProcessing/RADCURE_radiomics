library(MultiAssayExperiment)
library(SummarizedExperiment)
library(data.table)

getAssay <- function(experimentSEO, assayName) {
    assayN <- assay(experimentSEO, assayN)
    assayN_transposed <- t(assayN)
    return(assayN_transposed)
}

makeExperimentDataframe <- function(experimentSEO, assayNameList = NULL) {

    # If assayNameList is NULL, get all assays from experiment
    if (is.null(assayNameList)) {
        assayNameList <- names(assays(experimentSEO))
    }

    assayList <- assays(experimentSEO)
    # Initialize dataframe to append each assay to
    completeExpDataframe <- data.frame()

    for (assayName in assayNameList) {
        # Get data frame for this assay
        assayDataframe <- as.data.frame(assayList[assayName])

        # Get out shape features
        shapeFeatures <- assayDataframe[rownames(assayDataframe) %like% "shape_", ]
        shapeFeatures <- subset(shapeFeatures, select = -c(group, group_name))

        # If there are shape features present and the complete dataframe is empty,
        # add them to the front of the dataframe
        if (nrow(shapeFeatures) != 0) {
            if (nrow(completeExpDataframe) == 0) {
                completeExpDataframe <- rbind(completeExpDataframe, shapeFeatures)
            }

            # Remove shape features from assayDataframe
            assayDataframe <- assayDataframe[!(row.names(assayDataframe) %in% row.names(shapeFeatures)), ]
        }

        # Combine group_name (assay name) with features
        rownames(assayDataframe) <- paste(assayDataframe$group_name, rownames(assayDataframe), sep = "")

        # Remove the group and groupname columns
        assayDataframe <- subset(assayDataframe, select = -c(group, group_name))

        # Merge it to the main dataframe by the patient ID
        completeExpDataframe <- rbind(completeExpDataframe, assayDataframe)
    }
    # Transpose entire dataframe so columns are features and rows are patient IDs
    t_completeExpDataframe <- transpose(completeExpDataframe)
    colnames(t_completeExpDataframe) <- rownames(completeExpDataframe)
    rownames(t_completeExpDataframe) <- colnames(completeExpDataframe)

    # Add the metadata back to the front of the csv
    patientMeta <- colData(experimentSEO)
    dfpatientMeta <- as.data.frame(patientMeta)
    finalDf <- cbind(dfpatientMeta, t_completeExpDataframe)

    return(finalDf)
}

# MAE file path
maeFilePath <- "data/RADCURE_complete_radiomic_MAE_v2.rds"
# load MAE .rds file
radiomicMAE <- readRDS(maeFilePath)
# extract output dir name from MAE location to save csvs in the same place
outDir <- basename(dirname(maeFilePath))
# study name for output file creation
studyName <- paste(outDir, "RADCURE_complete", sep = "/")
# Get list of Experiment names
experimentNameList <- names(experiments(radiomicMAE))
# for each experiment
for (experimentName in experimentNameList) {
    # get out the summarized experiment object
    experimentSEO <- experiments(radiomicMAE)[[experimentName]]
    # combine the assay into one big dataframe
    experimentDF <- makeExperimentDataframe(experimentSEO)

    # save out as csv or xlsx
    outputFileName <- paste(studyName, experimentName, "features.csv", sep = "_")
    write.csv(experimentDF, outputFileName)
}
