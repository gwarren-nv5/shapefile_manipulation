@echo off

call "C:\ProgramData\Anaconda3\Scripts\activate.bat" C:\ProgramData\Anaconda3
call activate C:\Scripts\conda_env\pdxbridge

set in_shape=%~1

:: prompt the user for the remaining inputs
set /P out_crs="Enter EPSG code for output coordinate system: "
set /P buffer_distance="Enter buffer distance(to be buffered from outermost wire): "
set /P buffer_units="Enter buffer units (m/ft): "

set python_script_path="Z:\RESOURCES\Production\Bridge\00_Scripts\ZZZ_GW_working\shapefile_manipulation\shapfile_manipulations.py"

python %python_script_path% %in_shape% %out_crs% %buffer_distance% %buffer_units%
pause