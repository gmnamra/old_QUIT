#!/usr/bin/python

import Tkinter as Tk
import tkFileDialog, tkMessageBox
import subprocess
import os, glob, errno

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

class StatusBar(Tk.Frame):
    def __init__(self, master):
		Tk.Frame.__init__(self, master)
		self.label = Tk.Label(self, bd=1, relief=Tk.SUNKEN, anchor=Tk.W, width = 64)
		self.label.pack(fill=Tk.X)
		self.label.pack_propagate(False)
	
    def set(self, format, *args):
		self.label.config(text=format % args)
		self.label.update_idletasks()

    def clear(self):
		self.label.config(text="")
		self.label.update_idletasks()

class App:
	def __init__(self, master):
		frame = Tk.Frame(master)
		frame.pack()
		
		input = Tk.Frame(frame)
		input.grid(row = 0, column = 0)
		Tk.Label(input, text = 'Input:').grid(row = 0, sticky=Tk.E)
		Tk.Label(input, text = 'Output:').grid(row = 1, sticky=Tk.E)
		self.in_entry  = Tk.Entry(input, width = 64)
		self.out_entry = Tk.Entry(input, width = 64)
		self.in_entry.grid(row = 0, column = 1)
		self.out_entry.grid(row = 1, column = 1)
		self.in_button = Tk.Button(input, text = "...", command = self.find_in)
		self.in_button.grid(row = 0, column = 2)
		self.out_button = Tk.Button(input, text = "...", command = self.find_out)
		self.out_button.grid(row = 1, column = 2)
		
		options = Tk.Frame(frame)
		options.grid(row = 1, column = 0)
		self.study = Tk.IntVar()
		self.study.set(1)
		Tk.Checkbutton(options, text = "Whole study", variable = self.study).grid(row = 0, column = 0, sticky = Tk.W)
		self.spm_scale = Tk.IntVar()
		self.spm_scale.set(1)
		Tk.Checkbutton(options, text = "Scale for SPM", variable = self.spm_scale).grid(row = 0, column = 1, sticky = Tk.W)
		self.embed_procpar = Tk.IntVar()
		self.embed_procpar.set(1)
		Tk.Checkbutton(options, text = "Embed procpar", variable = self.embed_procpar).grid(row = 0, column = 2, sticky = Tk.W)
		self.gz = Tk.IntVar()
		self.gz.set(0)
		Tk.Checkbutton(options, text = "Compress", variable = self.gz).grid(row = 0, column = 3, sticky = Tk.W)
		self.ignore_scout = Tk.IntVar()
		self.ignore_scout.set(1)
		Tk.Checkbutton(options, text = "Ignore scouts", variable = self.ignore_scout).grid(row = 0, column = 4, sticky = Tk.W)
		self.echo_mode = Tk.IntVar()
		self.echo_mode.set(-1)
		Tk.Label(options, text = "Multi-echo:").grid(row = 1, column = 0, sticky = Tk.E)
		Tk.Radiobutton(options, text = "All", variable = self.echo_mode, value = -1).grid(row = 1, column = 1, sticky = Tk.W)
		Tk.Radiobutton(options, text = "Sum", variable = self.echo_mode, value = -2).grid(row = 1, column = 2, sticky = Tk.W)
		Tk.Radiobutton(options, text = "Select:", variable = self.echo_mode, value=0).grid(row = 1, column = 3, sticky = Tk.W)
		self.echo_select = Tk.Spinbox(options, from_ = 0, to=999999, width = 8) # Hopefully that is high enough!
		self.echo_select.grid(row = 1, column = 4, sticky = Tk.W)
		
		self.go_button = Tk.Button(options, text = "Convert", command = self.go)
		self.go_button.grid(row = 0, column = 5, rowspan = 2, columnspan = 2, sticky=Tk.N+Tk.S+Tk.E+Tk.W)
		self.go_text = StatusBar(frame)
		self.go_text.set("Ready")
		self.go_text.grid(row = 2, sticky=Tk.W + Tk.E)
		self.go_text.grid_propagate(False)
		
		menu = Tk.Menu(root)
		root.config(menu=menu)
		filemenu = Tk.Menu(menu)
		menu.add_cascade(label="File", menu=filemenu)
		filemenu.add_command(label="Exit", command=frame.quit)
		master.title("FDF to Nifti Conversion Tool")
		self.master = master
		
	def find_in(self):
		self.in_entry.delete(0, Tk.END)
		self.in_entry.insert(0, tkFileDialog.askdirectory(initialdir = "/data/blinded/OSIRIS/", mustexist = True))
	
	def find_out(self):
		self.out_entry.delete(0, Tk.END)
		self.out_entry.insert(0, tkFileDialog.askdirectory(initialdir = "~", mustexist = True))

	def go(self):
		(inpath, inext) = os.path.splitext(os.path.normpath(self.in_entry.get()))
		(indir, inbase) = os.path.split(os.path.normpath(inpath))
		outpath = self.out_entry.get()
		
		command = 'fdf2nii -v'
		if self.spm_scale.get():
			command = command + '-s 10.0 '
		
		if self.embed_procpar.get():
			command = command + '-p '
		
		if self.gz.get():
			command = command + '-z '
		
		if self.echo_mode.get() >= 0:
			command = command + '-e ' + self.echo_select.get() + ' '
		else:
			command = command + '-e ' + str(self.echo_mode.get()) + ' '
		
		if self.study.get():
			if inext == ".img":
				tkMessageBox.showwarning("Wrong folder", "You must select the parent folder when converting a whole study.")
				return
			outpath = outpath + '/' + inbase
			mkdir_p(outpath)
			command = command + '-o ' + outpath + '/ '
			
			scans = glob.glob(inpath + '/*.img')
			for s in scans:
				if (self.ignore_scout.get() and (os.path.basename(s).lower().find("scout") != -1)):
					pass
				else:
					command = command + s + ' '
		else:
			if inext != ".img":
				tkMessageBox.showwarning("Wrong Extension", "You must select the .img folder when converting a single image.")
				return
			command = command + '-o ' + outpath + '/ ' + inpath + inext
		
		#print command
		self.go_text.set("Starting...")
		self.master.config(cursor = "watch")
		self.master.update_idletasks()
		p = subprocess.Popen(command,
							 shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		for line in iter(p.stdout.readline, ''):
			self.go_text.set(line.rstrip())
			self.master.update_idletasks()
		self.master.config(cursor = "")
		self.go_text.set("Finished")

root = Tk.Tk()
app = App(root)
root.resizable(0, 0)
root.mainloop()
