from setuptools import setup

APP = ['MainWindow.py']
DATA_FILES = []
OPTIONS = {
    'plist': {'CFBundleShortVersionString':'0.1.0',},
    "py2app":
    {'iconfile':'images/Logo.icns',"argv_emulation": True, "includes": ["sip", "PyQt4._qt"], 'plist': {'CFBundleIconFile': 'Logo.icns'}}
}

setup(
  app= APP,
  name='MLSTEZ',
  data_files=DATA_FILES,
  options=OPTIONS,
  setup_requires=["py2app"]) 
