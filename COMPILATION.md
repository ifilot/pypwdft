# Compilation details

## Compilation and testing under Linux Debian / Ubuntu

Compile locally

```python
python3 -m build
```

and install it locally

```python
pip3 install dist/pypwdft-<version>-py3-none-any.whl
```

and finally test it

```python
python3 -m pytest tests/*
```

## Outputting content of Python Wheel file

```bash
unzip -l dist/*.whl
```