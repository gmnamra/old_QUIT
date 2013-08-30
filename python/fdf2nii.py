#!/usr/bin/python

import Tkinter as Tk
import tkFileDialog 
import subprocess

root = Tk.Tk()

class App:
	def __init__(self, master):
		frame = Tk.Frame(master)
		frame.pack()
		
		Tk.Label(frame, text = 'Input Folder:').grid(row = 0, sticky=Tk.E)
		Tk.Label(frame, text = 'Output Folder:').grid(row = 1, sticky=Tk.E)
		self.in_entry  = Tk.Entry(frame)
		self.out_entry = Tk.Entry(frame)
		self.in_entry.grid(row = 0, column = 1)
		self.out_entry.grid(row = 1, column = 1)
		
		self.in_button = Tk.Button(frame, text = "...", command = self.find_in)
		self.in_button.grid(row = 0, column = 2)
		self.out_button = Tk.Button(frame, text = "...", command = self.find_out)
		self.out_button.grid(row = 1, column = 2)

		self.go_button = Tk.Button(frame, text = "Convert", command = self.go)
		self.go_button.grid(row = 2)
		self.go_text = Tk.StringVar()
		Tk.Label(frame, textvariable = self.go_text).grid(row = 2, column = 1, columnspan = 2)
		
	def find_in(self):
		self.in_entry.delete(0, Tk.END)
		self.in_entry.insert(0, tkFileDialog.askdirectory())
	
	def find_out(self):
		self.out_entry.delete(0, Tk.END)
		self.out_entry.insert(0, tkFileDialog.askdirectory())

	def go(self):
		inname  = self.in_entry.get()
		outname = self.out_entry.get()
		self.go_text.set("Starting...")
		root.config(cursor = "wait")
		root.update()
		p = subprocess.Popen('fdf2nii -o ' + outname + '/ ' + inname,
							 shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		for line in iter(p.stdout.readline, ''):
			self.go_text.set(line)
			root.update()
		root.config(cursor = "")
		self.go_text.set("Finished")

app = App(root)
root.mainloop()