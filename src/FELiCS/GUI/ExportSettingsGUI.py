
class ExportSettingsGUI():
    def __init__(self,mainGUI):
        '''Initializing the Export GUI
        \t Input:
        \t -mainGUI: The FELiCS main GUI object
        '''
        import tkinter as tk
        from FELiCS.GUI.GUISettings import labelFont
        from tkinter import scrolledtext
        self.window =tk.Toplevel(mainGUI.window)
        self.Export=mainGUI.param.Export
        self.window.title("Numerics Settings")

        self.workDir=mainGUI.workDir
        # Mesh File
        ExportFolderColumn=1
        ExportFolderRow=1
        self.ExportFolderB = tk.Button(self.window, text='Choose export folder', command=self.chooseExportFolder)
        self.ExportFolderB.grid         (row = ExportFolderRow , column = ExportFolderColumn, rowspan = 1, columnspan = 1)
        self.ExportFolderStr    = tk.StringVar()
        self.ExportFolderL   = tk.Label(self.window, textvariable=self.ExportFolderStr)
        self.ExportFolderL.grid         (row = ExportFolderRow , column = ExportFolderColumn+1, rowspan = 1, columnspan = 1)
        self.ExportFolderStr.set(self.Export.ExportFolder)

        # Define export frame
        ExportFrameColumn=1
        ExportFrameRow=2
        self.ExportFrame = tk.LabelFrame(self.window, text="Export as...",font=labelFont())
        self.ExportFrame.grid        (row = ExportFrameRow, column = ExportFrameColumn, rowspan = 1, columnspan = 1, sticky="")

        # # Define vtk checkbox
        # self.VtkL = tk.Label(self.ExportFrame, text='...vtk', font=labelFont())
        # self.VtkL.grid(row = 1 , column = 1, rowspan = 1, columnspan = 1)
        # self.vtkBool=tk.BooleanVar()
        # self.vtkCheckBox=tk.Checkbutton(self.ExportFrame,
        #   text='',
        #   variable = self.vtkBool,
        #   command=self.doNothing)
        # self.vtkCheckBox.grid(row = 1 , column = 2, rowspan = 1, columnspan = 1)
        # self.vtkBool.set(self.Export.vtk)

        # # Define mat checkbox
        # self.matL = tk.Label(self.ExportFrame, text='...mat', font=labelFont())
        # self.matL.grid(row = 2 , column = 1, rowspan = 1, columnspan = 1)
        # self.matBool=tk.BooleanVar()
        # self.matCheckBox=tk.Checkbutton(self.ExportFrame,
        #   text='',
        #   variable = self.matBool,
        #   command=self.doNothing)
        # self.matCheckBox.grid(row = 2 , column = 2, rowspan = 1, columnspan = 1)
        # self.matBool.set(self.Export.mat)
        #
        # # Define hdf checkbox
        # self.hdfL = tk.Label(self.ExportFrame, text='...hdf', font=labelFont())
        # self.hdfL.grid(row = 3 , column = 1, rowspan = 1, columnspan = 1)
        # self.hdfBool=tk.BooleanVar()
        # self.hdfCheckBox=tk.Checkbutton(self.ExportFrame,
        #   text='',
        #   variable = self.hdfBool,
        #   command=self.doNothing)
        # self.hdfCheckBox.grid(row = 3 , column = 2, rowspan = 1, columnspan = 1)
        # self.hdfBool.set(self.Export.hdf)

        # Define Video checkbox
        self.VideoL = tk.Label(self.ExportFrame, text='...Video', font=labelFont())
        self.VideoL.grid(row = 4 , column = 1, rowspan = 1, columnspan = 1)
        self.VideoBool=tk.BooleanVar()
        self.VideoCheckBox=tk.Checkbutton(self.ExportFrame,
            text='',
            variable = self.VideoBool,
            command=self.doNothing)
        self.VideoCheckBox.grid(row = 4 , column = 2, rowspan = 1, columnspan = 1)
        self.VideoBool.set(self.Export.Video)

        # # Define h5 checkbox
        # self.h5L = tk.Label(self.ExportFrame, text='...h5', font=labelFont())
        # self.h5L.grid(row = 5 , column = 1, rowspan = 1, columnspan = 1)
        # self.h5Bool=tk.BooleanVar()
        # self.h5CheckBox=tk.Checkbutton(self.ExportFrame,
        #   text='',
        #   variable = self.h5Bool,
        #   command=self.doNothing)
        # self.h5CheckBox.grid(row = 5 , column = 2, rowspan = 1, columnspan = 1)
        # self.h5Bool.set(self.Export.h5)

        # Buttons to close and save and close
        self.CancelB     = tk.Button(self.window,text='Cancel', width=10, command=lambda: self.Cancel())
        self.CancelB.grid        (row = 99, column = 1, rowspan = 1, columnspan = 2, sticky="")
        self.SaveNCloseB     = tk.Button(self.window,text='Save&Close', width=10, command=lambda: self.SaveNClose(mainGUI))
        self.SaveNCloseB.grid        (row = 99, column = 3, rowspan = 1, columnspan = 2, sticky="")
        self.refresh()

    def ReturnParametersToMain(self,mainGUI):
        '''Passing the settings to the parameters object of the mainGUI
        \t Input:
        \t -mainGUI: Needed to pass the changes in the parameter file'''
        mainGUI.param.Export.ExportFolder=self.ExportFolderStr.get()
        # mainGUI.param.Export.vtk=self.vtkBool.get()
        # mainGUI.param.Export.mat=self.matBool.get()
        # mainGUI.param.Export.hdf=self.hdfBool.get()
        mainGUI.param.Export.Video=self.VideoBool.get()

    def chooseExportFolder(self):
        from os import path,makedirs
        from tkinter import filedialog
        # Get output directory from user. os.path.relpath provides the relative path to the file
        self.SolutionDirectory = path.relpath(filedialog.askdirectory(initialdir='./'),self.workDir)
        # Write the output directory to the GUI
        self.ExportFolderStr.set(self.SolutionDirectory)
        # If the directory does not exist, create it
        if not path.exists(self.SolutionDirectory):
            os.makedirs(self.SolutionDirectory)

    def Cancel(self):
        '''Function closing the current window'''
        self.window.destroy()

    def SaveNClose(self,mainGUI):
        ''' Function when user pushes Save and Close
        \tInput:
        \t-mainGUI: mainGUI object needed to refresh'''
        self.ReturnParametersToMain(mainGUI)
        self.window.destroy()
        mainGUI.refresh()

    def refresh(self,*args):
        '''Function refreshing the GUI, depending on the current inputs'''
        pass

    def doNothing(self,*args):
        pass
