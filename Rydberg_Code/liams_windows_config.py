### this file determines the location of files in my local system
class liams_windows_config(): 

    def file_system_params(self):
        self.project_loc = "C:/Users/liamw/Desktop/RYDBERG_GIT-1/src/PhD_Rydberg_code_version_1/Springer_code"
        
        pass

    def __str__(self):
        return("Project Loc: " + self.project_loc)

    def __init__(self): 
        self.file_system_params()
        pass
    
