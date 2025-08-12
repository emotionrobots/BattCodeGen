# BattCodeGen
Battery Model C Code Generator

## Installation
```
python3 -m venv myvenv 

source myvenv/bin/activate

cd myvenv

git clone https://github.com/emotionrobots/BattCodeGen.git

cd BattCodeGen
```


## Install Prerequisits
```
sudo apt-get update

sudo apt-get install build-essential libsundials-dev
```


## Generate the C code
```
python batt_codegen.py
```

You should see `batt_model.c` and `batt_mode.h` in the directory.


## Build the C code
```
mv batt_model.c src/.

mv batt_model.h include/.

make

```

## Run the C model
```
bin/run_model

```

