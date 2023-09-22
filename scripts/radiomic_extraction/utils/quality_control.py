import pydicom
import numpy as np
import os
import pandas
"""
QualityControl class is used to instantiate an object to run the following QC checks on your DICOM data:
Patient ID, axial resolution, slice thickness, modality, pixel spacing, number of scans, kernel and missing and overlapping slices. 
This code is based on QC checks done in precision-medicine-toolbox (https://github.com/primakov/precision-medicine-toolbox)
"""
class QualityControl:
    """Init QualityControl object.

    Parameters
    ----------
    qc_path
        Dictionary of type: qc_params = {
            'modality': '',
            'patient_ID': '',
            'kernel_list': list(),
            'slice_thickness_range': list(),
            'scan_length_range': list(),
            'spacing_range': list()
        }

    img_path
        Path to directory containing the DICOM series to load.

    Returns
    -------
    None. 

    """
    def __init__(self, qc_params=dict(), imgPath=str("")) -> None:
        self.scans = list()
        self.qc_params = qc_params
        self.scans = self.__read_scans(imgPath)
        
    """Reads all the individual scans in given directory. 

    Parameters
    ----------
    img_path
        Path to directory containing the DICOM series to load.

    Returns
    -------
    List of scans loaded via pydicom. 

    """
    def __read_scans(self, path):
        scans = []
        for s in os.listdir(path):
            temp_file = pydicom.read_file(os.path.join(path, s), force=True)
            scans.append(temp_file)
        return scans
    
    """Checks patient ID of all scans with given value. 

    Parameters
    ----------
    ID
        String of expected PatientID for each scan.

    Returns
    -------
    True or False if the check passed.  

    """
    def __check_patient_id(self, ID):
        for slice in self.scans:
            if slice.PatientID != ID:
                return False
        return True
    
    """Determines if all scans are axial resolution. 

    Parameters
    ----------
    None

    Returns
    -------
    True or False if the check passed.  

    """
    def __axial_resolution(self):
        for slice in self.scans:
            if not int(slice.ImageOrientationPatient[0]) and not int(slice.ImageOrientationPatient[4]) and int(
                slice.ImageOrientationPatient[1]) and int(slice.ImageOrientationPatient[2]) and not int(
                slice.ImageOrientationPatient[3]) and int(slice.ImageOrientationPatient[5]):
                return False
        return True
    
    """Checks if all scans have a SliceThickness in a certain range. 

    Parameters
    ----------
    thickness_range
        List of two elements indicating the min and max for slice thickness. 

    Returns
    -------
    True or False if the check passed.  

    """
    def __check_slice_thickness(self, thickness_range):
        res = True
        slice_thickness = [np.round(x.SliceThickness, 1) for x in self.scans]
        for s in slice_thickness:
            if not (s >= thickness_range[0] and s <= thickness_range[1]):
                res = False
        return res

    """Checks if all scans are the same modality. 

    Parameters
    ----------
    modality
        String for expected modality of each scan. 

    Returns
    -------
    True or False if the check passed.  

    """
    def __check_modality(self, modality):
        for slice in self.scans:
            if slice.Modality != modality:
                return False
        return True
    
    """Checks if all scans have a convolution kernel specified by the input. 

    Parameters
    ----------
    kernel
        List of expected kernels for each scan. 

    Returns
    -------
    True or False if the check passed.  

    """
    def __check_kernel(self, kernels=list()):
        kernels = [x.lower() for x in kernels]
        kernels_list = [x.ConvolutionKernel for x in self.scans]
        if kernels_list.count(kernels_list[0]) == len(kernels_list) and kernels_list[0].lower() in set(kernels):
            return True
        else:
            return False

    """Checks if any scans overlap or are missing.  

    Parameters
    ----------
    None

    Returns
    -------
    True or False if the check passed.  

    """
    def __check_missing_overlapping_slices(self):
        temp_spacing = []
        scan_real_range = np.round(list(float(x.ImagePositionPatient[2]) for x in self.scans), 1)
        for i in range(len(scan_real_range) - 1):
            temp_spacing.append(np.round(scan_real_range[i + 1] - scan_real_range[i], 1))

        if temp_spacing.count(temp_spacing[0]) == len(temp_spacing):
            return True
        else:
            return False
        
    """Checks if all scans have the specified pixel spacing.

    Parameters
    ----------
    spacing_range
        List of size two with min and max for PixelSpacing range. 

    Returns
    -------
    True or False if the check passed.  

    """
    def __check_pixel_spacing(self, spacing_range):
        for scan in self.scans:
            ps = scan.PixelSpacing
            if not (ps[0]==ps[1] and ps[0]>=spacing_range[0] and ps[1]<=spacing_range[1]):
                return False
        return True
    
    """Checks if all scans are in the length range specified by input. 

    Parameters
    ----------
    scan_length_range
        List of size two with min and max for the expected range for the number of scans. 

    Returns
    -------
    True or False if the check passed.  

    """
    def __check_scan_length_range(self, scan_length_range=list()):
        return len(self.scans)>=scan_length_range[0] and len(self.scans) <=scan_length_range[1]

    """Runs all the checks specified above. 

    Parameters
    ----------
    None

    Returns
    -------
    Boolean pandas dataframe determining the outcome of each check. 

    """
    def run_checks(self):
        modality_check = self.__check_modality(self.qc_params['modality'])
        patient_id_check = self.__check_patient_id(self.qc_params['patient_ID'])
        axial_resolution_check = self.__axial_resolution()
        missing_overlapping_slices_check = self.__check_missing_overlapping_slices()
        scan_length_range_check = self.__check_scan_length_range(self.qc_params['scan_length_range'])
        pixel_spacing_check = self.__check_pixel_spacing(self.qc_params['spacing_range'])
        slice_thickness_check = self.__check_slice_thickness(self.qc_params['slice_thickness_range'])
        kernel_check = self.__check_kernel(self.qc_params['kernel_list'])
        return pandas.DataFrame([modality_check, patient_id_check, axial_resolution_check, missing_overlapping_slices_check, scan_length_range_check, pixel_spacing_check, slice_thickness_check, kernel_check])
    