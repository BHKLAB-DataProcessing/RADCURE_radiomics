import os
import pandas as pd

def saveDataframeCSV(dataframe: pd.DataFrame,
                     outputFilePath: str):
    ''' Function to save a pandas Dataframe as a csv file with the index removed.
        Checks if the path in the outputFilePath exists and will create any missing directories.

    Parameters
    ----------
    dataframe
        Pandas dataframe to save out as a csv
    outputFilePath
        Full file path to save the dataframe out to.
    '''

    if not outputFilePath.endswith('.csv'):
        raise ValueError("This function saves .csv files, so outputFilePath must end in .csv")

    # Make directory if it doesn't exist
    if not os.path.exists(os.path.dirname(outputFilePath)):
        os.makedirs(os.path.dirname(outputFilePath))

    # Save out feature set
    dataframe.to_csv(outputFilePath, index=False)

    return