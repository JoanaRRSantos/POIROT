#import tkinter as tk
import customtkinter as ctk
import appDevelopmentNCBI 
#from tkinter import ttk
#from tkinter import filedialog
from customtkinter import filedialog
from CTkMessagebox import CTkMessagebox

class MyTabView(ctk.CTkTabview):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        self.folderDirectory = ""
        self.folderDirectoryToSaveFiles=""
        self.fileName=""
    
        # create tabs
        self.add("Add Input File")
        #self.add("Run")

        # add widgets on tabs
        #ADD INPUT FILE TAB
        #label1
        self.label1 = ctk.CTkLabel(master=self.tab("Add Input File"), text="Which directory you have the input file?")
        self.label1.grid(row=0, column=0, padx=20, pady=10)
        #button1 - button to submit the folder where the input file is
        self.button1 = ctk.CTkButton(master=self.tab("Add Input File"), text='Search', command= self.getdirctoryteste)
        self.button1.grid(row=0, column=1)
        # label2
        self.label2 = ctk.CTkLabel(master=self.tab("Add Input File"), text="Which directory you want to save the results?")
        self.label2.grid(row=1,column=0)

        # button2 - button to submit the folder where to save the results
        self.button2 = ctk.CTkButton(master=self.tab("Add Input File"), text='Search', command=self.getFolderSavingDirectory)
        self.button2.grid(row=1, column=1)
        
        # label3
        self.label3 = ctk.CTkLabel(master=self.tab("Add Input File"), text="Tell me which file is the input file?")
        self.label3.grid(row=2,column=0)
        
        # button3 - button to submit the file
        self.button3 = ctk.CTkButton(master=self.tab("Add Input File"), text='Search', command = self.getFile)
        self.button3.grid(row=2, column=1)
        
        #label4 - to give a space between buttons
        self.label4 = ctk.CTkLabel(master=self.tab("Add Input File"), text=" ")
        self.label4.grid(row=3,column=0)
        #button4 - to run the program that does the BLAST automaticly
        
        self.button4 = ctk.CTkButton(master=self.tab("Add Input File"), text='Run BLAST', command = self.runProgram)
        self.button4.grid(row=4, column=1)
        """
        progress bar do program run 
        self.progressbar = ctk.CTkProgressBar(master=self.tab("Run"))
        self.progressbar.grid(padx=20, pady=10)  
        """
    #functions for buttons:
    def getdirctoryteste(self):
        #folderDirectory=t.get()
        tempDirectory = filedialog.askdirectory(initialdir = ".", title = "Open File")
        self.folderDirectory =  r'%s' % tempDirectory
        print(self.folderDirectory)
        
    def getFolderSavingDirectory(self):
        #folderDirectory=t.get()
        tempfolderDirectoryToSaveFiles = filedialog.askdirectory(initialdir = ".", title = "Open File")
        self.folderDirectoryToSaveFiles = r'%s' % tempfolderDirectoryToSaveFiles
        print(self.folderDirectoryToSaveFiles)

    def getFile(self):
        #folderDirectory=t.get()
        tempfileName = filedialog.askopenfilename(initialdir = ".", title = "Open File")
        self.fileName = r'%s' % tempfileName
        print(self.fileName)
    
    def runProgram(self):
        teste = appDevelopmentNCBI.compareSequencesNCBI(directory= self.folderDirectory, folderDirectoryToSaveFiles= self.folderDirectoryToSaveFiles, fileName = self.fileName)
        teste.getDirectoryFromUser()
        teste.getFileNameFromUser()
        teste.processInputFile()
        teste.saveBLAST()
        teste.processBLASTresults()
        teste.constructResultsFile()
        #teste.removeIntermediateFiles()
        ctk.CTkButton(master=self.tab("Add Input File"), text="Program Done", command=self.show_checkmark).grid(padx=20, pady=10, sticky="news")
    
    def show_checkmark(self):
    # Show some positive message with the checkmark icon
        CTkMessagebox(title="Info", message="The results are ready in your folder! Thank you for using SequencesComparererer.", icon="check", option_1="Exit")



class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        # configure the root window
        self.title('SequencesComparererer')
        self.tab_view = MyTabView(master=self)
        self.tab_view.grid(row=0, column=0)
        """
        self.teste_menu = tk.Menubutton(self, text="Type of BLAST:")
        self.teste_menu.menu = tk.Menu(self.teste_menu, tearoff=0)
        self.teste_menu["menu"] = self.teste_menu.menu
        
        self.teste_menu.menu.add_command(label="blastn", command = lambda: print("This is the BLASTN option"))
        self.teste_menu.menu.add_command(label="blastp", command = lambda: print("This is the BLASTP option"))
        self.teste_menu.menu.add_command(label="blastx", command = lambda: print("This is the BLASTX option"))
        self.teste_menu.grid(row=4, column=0)
        
        #button4 - to run the program that does the BLAST automaticly
        self.button4 = ctk.CTkButton(self, text='Run BLAST', command = self.runProgram)
        self.button4.grid(row=5, column=2)
    
    """

        
if __name__ == "__main__":
    app = App()
    app.mainloop()