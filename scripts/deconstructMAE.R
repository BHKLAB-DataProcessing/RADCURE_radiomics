
# load MAE .rds file
radiomicMAE <- readRDS("data/RADCURE_complete_radiomic_MAE.rds")
# Get list of Experiment names
experimentNameList <- names(experiments(radiomicMAE))
# for each experiment
for (experimentName in experimentNameList) {
    # get out the summarized experiment object
    experimentSEO <- experiments(radiomicMAE)[[experimentName]]
    # combine the assay into one big dataframe

    # coldata at the front of the dataframe

    # save out as csv or xlsx

}
    