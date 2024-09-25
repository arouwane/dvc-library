echo " "
echo "############################### File setup_imageRoutines.py ###############################"
echo " " 
python setup_imageRoutines.py build_ext --inplace 
echo "############################### File setup_toolsRoutines.py ###########################"  
python setup_toolsRoutines.py build_ext --inplace  
echo "############################### File setup_assemblyRoutines.py ###########################"  
python setup_assemblyRoutines.py build_ext --inplace  
#echo "############################### File setup_assembly_Bspline_3d.py ###########################"  
#python setup_assembly_Bspline_3d.py build_ext --inplace  
