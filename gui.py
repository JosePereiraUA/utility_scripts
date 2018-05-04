#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Tkinter import *
import ttk
import tkFileDialog as filedialog
import tkMessageBox as alert
import os

class Global:

    LABEL_WIDTH     = 300
    FILE_WIDTH      = 50
    BACKGROUND      = "white"
    PAD             = 10
    FONT_SIZE       = 10
    FONT            = "System"
    HIGHLIGHT_COLOR = "steel blue"
    STICKY          = "w"
    BUTTON_SIZE     = 12
    RELIEF          = SUNKEN
    BORDER_WIDTH    = 1


class Input:

    DEFAULT_LABEL_WIDTH = Global.LABEL_WIDTH
    DEFAULT_PAD         = Global.PAD
    DEFAULT_WIDTH       = 15
    DEFAULT_FONT_SIZE   = Global.FONT_SIZE
    DEFAULT_STICKY      = Global.STICKY
    DEFAULT_BACKGROUND  = Global.BACKGROUND
    DEFAULT_FONT        = Global.FONT

    def __init__(self, master, label, default_value = "...", row = 0, column = 0, columnspan = 1,
        padx = DEFAULT_PAD, pady = DEFAULT_PAD, width = DEFAULT_WIDTH, font_size = DEFAULT_FONT_SIZE,
        font = DEFAULT_FONT, sticky = DEFAULT_STICKY, background = DEFAULT_BACKGROUND,
        label_width = DEFAULT_LABEL_WIDTH):
        
        self.value = StringVar()
        self.value.set(default_value)
        
        frame = Frame(master, background = background)
        frame.grid(row = row, column = column, columnspan = columnspan, sticky = sticky)
        frame.grid_columnconfigure(0, minsize = label_width)
        
        label = ttk.Label(frame, text = label)
        label.grid(row = row, column = 0, padx = padx, pady = pady, sticky = 'e')
        label.config(font=(font, font_size))

        ttk.Entry(frame, width = width, textvariable = self.value)\
            .grid(row = row, column = 1, padx = padx)


    def __call__(self):
        return self.value.get()


class CheckBox:

    DEFAULT_PAD         = Global.PAD
    DEFAULT_FONT_SIZE   = Global.FONT_SIZE
    DEFAULT_STICKY      = Global.STICKY
    DEFAULT_FONT        = Global.FONT
    DEFAULT_HOVER_COLOR = Global.HIGHLIGHT_COLOR
    DEFAULT_COLOR       = Global.BACKGROUND

    def __init__(self, master, label, default_status = True, row = 0, column = 0, 
        padx = DEFAULT_PAD, pady = DEFAULT_PAD, hover_color = DEFAULT_HOVER_COLOR,
        font_size = DEFAULT_FONT_SIZE, font = DEFAULT_FONT, sticky = DEFAULT_STICKY,
        color = DEFAULT_COLOR):
        
        self.status = IntVar(value = default_status)
        check_box = Checkbutton(master, text=' ' + label, variable=self.status,
            activebackground = hover_color, background = color, highlightthickness = 0, bd = 0)
        check_box.grid(row = row, column = column, padx = padx, pady = pady, sticky = sticky)
        check_box.config(font=(font, font_size))
        if default_status: check_box.select()


    def __call__(self):
        return self.status.get()


class Title:
    
    DEFAULT_PAD       = Global.PAD
    DEFAULT_FONT_SIZE = Global.FONT_SIZE
    DEFAULT_STICKY    = ''
    DEFAULT_FONT      = Global.FONT
    DEFAULT_COLOR     = Global.BACKGROUND


    def __init__(self, master, label, row = 0, padx = DEFAULT_PAD, pady = DEFAULT_PAD, columnspan = 1,
        font_size = DEFAULT_FONT_SIZE, font = DEFAULT_FONT, sticky = DEFAULT_STICKY, color = DEFAULT_COLOR):
        
        label = ttk.Label(master, text = label, background = color)
        label.grid(row = row, columnspan = columnspan, padx = padx, pady = pady, sticky = sticky)
        label.config(font=(font, font_size))


class SelectFile:

    DEFAULT_LABEL_WIDTH = Global.LABEL_WIDTH
    DEFAULT_PAD         = Global.PAD
    DEFAULT_WIDTH       = Global.FILE_WIDTH
    DEFAULT_FONT_SIZE   = Global.FONT_SIZE
    DEFAULT_BUTTON_SIZE = Global.BUTTON_SIZE
    DEFAULT_STICKY      = Global.STICKY
    DEFAULT_BACKGROUND  = Global.BACKGROUND
    DEFAULT_FONT        = Global.FONT
    DEFAULT_HOVER_COLOR = Global.HIGHLIGHT_COLOR

    def __init__(self, master, label, default_value = "Please select a file ...",
        row = 0, column = 0, padx = DEFAULT_PAD, pady = DEFAULT_PAD, width = DEFAULT_WIDTH,
        font_size = DEFAULT_FONT_SIZE, font = DEFAULT_FONT, sticky = DEFAULT_STICKY,
        button_size = DEFAULT_BUTTON_SIZE, hover_color = DEFAULT_HOVER_COLOR,
        background = DEFAULT_BACKGROUND, label_width = DEFAULT_LABEL_WIDTH, search_directory_only = False):
        
        self.value = StringVar()
        self.value.set(default_value)
        self.dir = search_directory_only
        
        frame = Frame(master, background = background)
        frame.grid(row = row, column = column, columnspan = 2, sticky = sticky)
        frame.grid_columnconfigure(0, minsize = label_width)

        label = ttk.Label(frame, text = label)
        label.grid(row = row, column = 0, padx = padx, pady = pady, sticky = 'e')
        label.config(font=(font, font_size))

        ttk.Entry(frame, width = width, textvariable = self.value)\
            .grid(row = row, column = 1, padx = padx)

        Button(frame, text="Browse", width = button_size, command= self.callback,
            activebackground = hover_color, relief = 'flat', background = background)\
            .grid(row = row, column = 2, padx = padx)


    def callback(self):
        if self.dir:
            self.value.set(filedialog.askdirectory())
        else:
            self.value.set(filedialog.askopenfilename())


    def __call__(self):
        return self.value.get()


class Container:

    DEFAULT_FONT_SIZE    = Global.FONT_SIZE
    DEFAULT_BACKGROUND   = Global.BACKGROUND
    DEFAULT_FONT         = Global.FONT
    DEFAULT_HOVER_COLOR  = Global.HIGHLIGHT_COLOR
    DEFAULT_RELIEF       = Global.RELIEF
    DEFAULT_BORDER_WIDTH = Global.BORDER_WIDTH
    DEFAULT_STICKY       = Global.STICKY


    def __init__(self, master, default_text = '', row = 0, column = 0, background = DEFAULT_BACKGROUND,
        font_color = DEFAULT_HOVER_COLOR, font = DEFAULT_FONT, hover_color = DEFAULT_HOVER_COLOR,
        relief = DEFAULT_RELIEF, borderwidth = DEFAULT_BORDER_WIDTH, font_size = DEFAULT_FONT_SIZE,
        sticky = DEFAULT_STICKY):

        self.frame = LabelFrame(master, text = default_text, foreground = font_color,
                background = background, borderwidth = borderwidth, relief = relief)
        self.frame.grid(row = row, column = column, padx = 10)
        self.frame.config(font=(font, font_size))

    def __call__(self):
        return self.frame


class Terminal:

    DEFAULT_PAD         = Global.PAD
    DEFAULT_WIDTH       = 800
    DEFAULT_BUTTON_SIZE = 12
    DEFAULT_FONT_SIZE   = Global.FONT_SIZE
    DEFAULT_BACKGROUND  = Global.BACKGROUND
    DEFAULT_FONT_COLOR  = Global.HIGHLIGHT_COLOR
    DEFAULT_FONT        = Global.FONT

    def __init__(self, master, default_text = '...', row = 0, background = DEFAULT_BACKGROUND,
        width = DEFAULT_WIDTH, font_size = DEFAULT_FONT_SIZE, font = DEFAULT_FONT,
        font_color = DEFAULT_FONT_COLOR, pady = DEFAULT_PAD, columnspan = 1):

        self.terminal = Message(master, fg = font_color, text = default_text, width = width,
            background = background)
        self.terminal.grid(row = row, pady = pady, columnspan = columnspan)
        self.terminal.config(font=(font, font_size))


    def __call__(self, message):
        self.terminal.config(text = message)


class Footer:

    DEFAULT_PAD         = Global.PAD
    DEFAULT_FONT_SIZE   = Global.FONT_SIZE
    DEFAULT_STICKY      = Global.STICKY
    DEFAULT_BACKGROUND  = Global.BACKGROUND
    DEFAULT_FONT        = Global.FONT
    DEFAULT_FONT_COLOR  = Global.HIGHLIGHT_COLOR

    def __init__(self, master, text, row = 0, columnspan = 1, sticky = DEFAULT_STICKY,
        font_color = DEFAULT_FONT_COLOR, background = DEFAULT_BACKGROUND, font = DEFAULT_FONT,
        font_size = DEFAULT_FONT_SIZE, pad = DEFAULT_PAD):

        footer = ttk.Label(master, text = text, foreground = font_color, background = background)
        footer.grid(row = row, columnspan = columnspan, sticky = 'se', pady = pad, padx = pad)
        footer.config(font=(font, font_size))


class Slice:

    DEFAULT_LABEL_WIDTH = Global.LABEL_WIDTH
    DEFAULT_WIDTH       = 7
    DEFAULT_BACKGROUND  = Global.BACKGROUND
    DEFAULT_STICKY      = Global.STICKY

    def __init__(self, master, text_from, default_from, text_to, default_to, row = 0, column = 0,
        columnspan = 1, width = DEFAULT_WIDTH, background = DEFAULT_BACKGROUND,
        label_width = DEFAULT_LABEL_WIDTH, sticky = DEFAULT_STICKY):

        frame = Frame(master, background = background)
        frame.grid(row = row, column = column, columnspan = columnspan, sticky = sticky)

        self.start = Input(frame, text_from, default_from, width = width, label_width = label_width)
        self.end   = Input(frame, text_to  , default_to, column = 1, width = width, label_width = 0)


class BasicGUIFunctions:

    def run(self, command, update_progress_bar = True):
        """
        Instantiate a new terminal command (GROMACS)
        """
        command = command + " " + self.debug
        os.popen(command)
        if update_progress_bar:
            self.update_progress_bar()


    def error(self, error_title, message, auto_exit = True):
        """
        Displays a pop-up alert to user and exits program.
        """

        alert.showerror(error_title, message)
        if auto_exit:
            exit(1)


    def update_progress_bar(self):
        """
        Increases the completness of the progress bar by one step.
        """

        self.progress_bar['value'] = self.progress_bar['value'] + 1
        self.root.update_idletasks()

    
    def update_terminal(self, message):
        """
        Shows a message in the GUI terminal.
        """

        self.terminal(message)
        self.root.update_idletasks()

    
    def clean(self):
        """
        Delete unnecessary files.
        """

        self.update_terminal("Cleaning unnecessary files ...")
        self.run("rm -rf *#")
        self.update_terminal("All tasks successefully performed.")

    
    def add_basic_utilities(self, master, min_row,
        pad = 10, background = 'white', hover_color = 'steel blue', footer = 'José Pereira @ 2018'):

        Button(master, text="Start MD",
            width = 50, command= self.process, activebackground = hover_color,
            relief = 'flat', background = background)\
            .grid(row = min_row, padx = pad, columnspan = 2)
        self.progress_bar = ttk.Progressbar(master, orient = 'horizontal', mode = 'determinate', length = 400)
        self.progress_bar.grid(columnspan = 2, row = min_row + 1, pady = pad)
        self.terminal     = Terminal(master, row = min_row + 2, columnspan = 2)
        self.footer       = Footer(master, footer, row = min_row + 3, columnspan = 2)

    
    def config_master(self, master,
        title = "Tool", version = "DEV_MODE", resizable = True, background = 'white'):

        master.configure(background = background)
        master.title("{title} v{version}".format(title = title, version = version))
        if not resizable: master.resizable(0, 0)
        self.root = master


class Step:

    def __init__(self, title):
        self.title = title
        self.mdp   = title + '.mdp'
        self.tpr   = title + '.tpr'
        self.n_steps     = 0
        self.print_every = 0

    def config_mdp(self, mdp_dir):

        lines = []
        
        with open(mdp_dir + '/' + self.mdp, 'r') as file_in:
            for line in file_in:
                if line.startswith("nsteps"):
                    line = line[:-2]
                    line = line + " " + str(self.n_steps) + '\n'
                elif line.startswith("nstxout") or line.startswith("nstenergy") or line.startswith("nstlog") or line.startswith("nstxout-compressed"):
                    line = line[:-2]
                    line = line + " " + str(self.print_every) + '\n'
                lines.append(line)
        
        with open(self.mdp, 'w') as file_out:
            for line in lines:
                file_out.write(line)
