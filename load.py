from tkinter import *

splash_root = Tk()
splash_root.title("COLEFIT")
splash_root.geometry("300x200")
splash_root.overrideredirect(True) #hide title bar

splash_label = Label(splash_root, text="test", font=("Helvetica", 18))
splash_label.pack(pady=20)

def main_window():
    splash_root.destroy()
    root = Tk()
    root.title("COLEFIT")
    root.geometry("500x300")

splash_root.after(3000, main_window)
mainloop()
