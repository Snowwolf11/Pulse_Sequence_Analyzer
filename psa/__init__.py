from .GUI import main as tkinter_main


def main(gui="tkinter"):
    if gui == "tkinter":
        return tkinter_main()

    if gui == "dash":
        from .GUI_dash import main as dash_main
        return dash_main()

    if gui == "pyside6":
        from .GUI_pyside6 import main as pyside6_main
        return pyside6_main()

    raise ValueError(f"Unknown GUI option: {gui}")
