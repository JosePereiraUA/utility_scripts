from Tkinter import *
import ttk
import tkFileDialog as filedialog

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