import argparse
import psa

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gui",
        choices=["tkinter", "dash", "pyside6"],
        default="pyside6",
    )
    args = parser.parse_args()

    psa.main(gui=args.gui)