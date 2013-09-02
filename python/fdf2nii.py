#!/usr/bin/python

import Tkinter as Tk
import tkFileDialog, tkMessageBox
import subprocess
import os

class App:
	def __init__(self, master):
		frame = Tk.Frame(master)
		frame.pack()
		
		input = Tk.Frame(frame)
		input.grid(row = 0, column = 0)
		Tk.Label(input, text = 'Input Folder:').grid(row = 0, sticky=Tk.E)
		Tk.Label(input, text = 'Output Folder:').grid(row = 1, sticky=Tk.E)
		self.in_entry  = Tk.Entry(input)
		self.out_entry = Tk.Entry(input)
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
		Tk.Checkbutton(options, text = "Whole study", variable = self.study).grid(row = 0, column = 0)
		self.spm_scale = Tk.IntVar()
		self.spm_scale.set(1)
		Tk.Checkbutton(options, text = "Scale for SPM", variable = self.spm_scale).grid(row = 0, column = 1)
		self.embed_procpar = Tk.IntVar()
		self.embed_procpar.set(1)
		Tk.Checkbutton(options, text = "Embed procpar", variable = self.embed_procpar).grid(row = 1, column = 0)
		self.gz = Tk.IntVar()
		self.gz.set(0)
		Tk.Checkbutton(options, text = "Compress", variable = self.gz).grid(row = 1, column = 1)
		
		go = Tk.Frame(frame)
		go.grid(row = 2)
		self.go_button = Tk.Button(go, text = "Convert", command = self.go)
		self.go_button.grid(row = 0, column = 2)
		self.go_text = Tk.StringVar()
		self.go_text.set("Ready")
		Tk.Label(go, textvariable = self.go_text, width=25).grid(row = 0, column = 0, columnspan = 2)
		
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
		self.out_entry.insert(0, tkFileDialog.askdirectory(initialdir = "~", mustexist = False))

	def go(self):
		inpath  = self.in_entry.get()
		outpath = self.out_entry.get()
		command = 'fdf2nii -o ' + outpath + '/ '
		if self.spm_scale.get():
			command = command + '-s 10.0 '
		if self.embed_procpar.get():
			command = command + '-p '
		if self.gz.get():
			command = command + '-z '
		if self.study.get():
			if inpath.endswith(".img"):
				tkMessageBox.showwarning("Wrong folder",
				                         "You must select the parent folder when converting a whole study.")
				return
			command = command + inpath + '/*.img'
		else:
			if not(inpath.endswith(".img")):
				tkMessageBox.showwarning("Wrong Extension",
										 "You must select the .img folder when converting a single image.")
				return
			command = command + inpath
		self.go_text.set("Starting...")
		self.master.config(cursor = "watch")
		self.master.update()
		p = subprocess.Popen(command,
							 shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		for line in iter(p.stdout.readline, ''):
			self.go_text.set(line)
			print line,
			self.master.update()
		self.master.config(cursor = "")
		self.go_text.set("Finished")

root = Tk.Tk()
app = App(root)
root.mainloop()
