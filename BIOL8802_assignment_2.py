#Assignment 2
import subprocess
import os

#Create a new class called FileManager
class FileManager:
#This class will assist you in downloading and uploading data from the cloud onto the local computer you are using.

#This class will need to be able to:
#(1) Take a local filename and convert it to a cloud filename
#   e.g., '/Users/pmcgrath7/TempData/Python3Class/test.txt' -> 'ptm_dropbox:/COS/BioSci/BioSci-McGrath/Python3Class/test.txt'
#(2) Download data from the cloud onto your computer
#(3) Upload data from your computer onto the cloud

#The class should have the following methods:

#__init__ -> This should store the name of the remote, the name of the master local directory, and the name of the master cloud directory
    def __init__(self, remote_name, mstr_local_dir, mstr_cld_dir):
        self.remote_name = remote_name
        self.mstr_local_dir = mstr_local_dir
        self.mstr_cld_dir = mstr_cld_dir

#convertCloudToLocal(filename) -> Should convert the name of the cloud location (which will include the remote) to the local name
    def convertCloudToLocal(self, filename):
        #help from https://www.geeksforgeeks.org/python-program-to-get-the-file-name-from-the-file-path/
        endname = os.path.basename(filename)
        local_name = mstr_local_dir + endname
        print('your new file name is' + local_name)
        return local_name

#convertLocalToCloud(filename) -> Should convert the name of the local location to the cloud location
    def convertLocalToCloud(self, filename):
        endname = os.path.basename(filename)
        cld_name = mstr_cld_dir + endname
        print('your new file name is' + cld_name)
        return cld_name

#uploadData(filename) -> uploads a local file to the cloud
    def uploadData(self, filename):
        subprocess.run(["rclone", "copy", local_name, cld_name], capture_output=True)
        print(filename + ' uploaded to ' + mstr_cld_dir)

#downloadData(filename) -> downloads a cloud file to the local drive
    def downloadData(self, filename):
        subprocess.run(["rclone", "copy", cld_name, local_name], capture_output=True)
        print (filename +' downloaded to ' + mstr_local_dir)