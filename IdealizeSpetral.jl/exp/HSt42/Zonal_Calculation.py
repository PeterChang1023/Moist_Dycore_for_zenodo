# import numpy as np
# import h5py

# class Zonal_Calculation:
#     def __init__(self, zonal_dimension: int, data_name: str):
#         self.idx = zonal_dimension
#         self.data_name = data_name

#     def zonal_mean(self, data):
#         return data.mean(axis=self.idx)

#     def zonal_anomaly(self, data):
#         mean = self.zonal_mean(data)
#         np.subtract(data, mean[:,:,:,np.newaxis], out=data)  # In-place subtraction
#         return data

#     def load_data(self, file_path):
#         with h5py.File(file_path, "r") as file:
#             # Directly return the numpy array
#             return np.array(file[self.data_name], dtype=np.float32)
# ChatGPT
import numpy as np
import h5py

class Zonal_Calculation:
    def __init__(self, zonal_dimension: int, data_name: str):
        self.idx = zonal_dimension
        self.data_name = data_name

    def zonal_mean(self, data):
        return np.mean(data, axis=self.idx, keepdims=True)

    def time_mean(self, data):
        return np.mean(data, axis=0, keepdims=True)

    def mean_exclude_latitude(self, data):
        return np.nanmean(data, axis=(0, 1, 3))

    def mean_exclude_latitude_for_prec(self, data):
        return np.sum(np.mean(data, axis=(1, 3)))

    def zonal_anomaly(self, data):
        mean = self.zonal_mean(data)
        np.subtract(data, mean, out=data)  # In-place subtraction
        return data

    def load_data(self, file_path):
        with h5py.File(file_path, "r") as file:
            return np.array(file[self.data_name], dtype=np.float32)