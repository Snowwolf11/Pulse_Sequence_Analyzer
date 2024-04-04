## Run from python
```
# Create a virtual python envrionment to not polute the global one.
# If you are already in a virtual environment or want this in your
# global one, just skip this
python -m venv .env
source .env/bin/activate

# Install the tool that you'll use to build the package
pip install maturin

# Build the python package
maturin build --release --features python
# OR build it with paralization enabled
maturin build --release --features python,par

# The generated whl file can then be installed using `pip`. Just google it :)
```

## Run from CLI
```
# Get help (The output will tell you how to provide the args)
# 
cargo run --release --features cli -- --help

#
```
