import os

folder_path = '/home/vfornemo/project/SPG_analyzer/tests/testset/'

for filename in os.listdir(folder_path):
    if filename.endswith('.mol'):
        new_filename = filename.split('_', 1)[-1]
        os.rename(os.path.join(folder_path, filename), os.path.join(folder_path, new_filename))