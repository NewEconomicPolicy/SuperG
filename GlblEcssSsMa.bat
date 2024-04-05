rem
@set setup_dir=E:\AbUniv\GlobalEcosseSuite\setup\
@set source_dir=G:\AbUnivGit\SuperGSpGlEc\GlblEcosseSiteSpecMa\
@set envmdlng_dir=G:\AbUnivGit\EnvMdllngModuls\
@set python_exe=E:\Python38\python.exe
@set PYTHONPATH=%envmdlng_dir%EnvModelModules;%envmdlng_dir%GlblEcosseMisc

@set initial_working_dir=%cd%
@chdir /D %setup_dir%site_spec
start cmd.exe /k "%python_exe% -W ignore %source_dir%GlblEcsseHwsdGUI.py"
@chdir /D %initial_working_dir%
