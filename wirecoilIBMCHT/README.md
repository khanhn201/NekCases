


# Dependencies (from Aya Hegazy)
## Embree
```
cd $HOME/local/src
git clone https://urldefense.com/v3/__https://github.com/RenderKit/embree.git__;!!DZ3fjg!74F3VvkmQkbCPdS1x-FPxcyi7MKN1RdLnNHY-C5KGJeNCjcyzw8Yn2NlJLgENMXHyNUiyA82qzTOTc1TSRUsMMm8Q1X8$ 
cd embree
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/local -DEMBREE_ISPC_SUPPORT=OFF -DEMBREE_TUTORIALS=OFF -DEMBREE_TASKING_SYSTEM=INTERNAL
cmake --build . -j 8
cmake --install .
```
## VTK
```
cd $HOME/local/src
git clone https://urldefense.com/v3/__https://github.com/Kitware/VTK.git__;!!DZ3fjg!74F3VvkmQkbCPdS1x-FPxcyi7MKN1RdLnNHY-C5KGJeNCjcyzw8Yn2NlJLgENMXHyNUiyA82qzTOTc1TSRUsMNDRtr3v$ 
cd VTK
git checkout v9.3.0
mkdir -p build
cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$HOME/local \
  -DVTK_BUILD_TESTING=OFF \
  -DVTK_BUILD_EXAMPLES=OFF \
  -DVTK_WRAP_PYTHON=OFF \
  -DVTK_WRAP_JAVA=OFF \
  -DVTK_GROUP_ENABLE_Rendering=WANT \
  -DVTK_GROUP_ENABLE_StandAlone=WANT \
  -DVTK_MODULE_ENABLE_VTK_IOGeometry=YES \
  -DVTK_MODULE_ENABLE_VTK_IOXML=YES \
  -DVTK_MODULE_ENABLE_VTK_IOLegacy=YES \
  -DVTK_MODULE_ENABLE_VTK_FiltersCore=YES \
  -DVTK_MODULE_ENABLE_VTK_FiltersGeneral=YES \
  -DVTK_MODULE_ENABLE_VTK_FiltersHybrid=YES \
  -DVTK_MODULE_ENABLE_VTK_RenderingCore=YES

cmake --build . -j 8
cmake --install .
```
There is a small bug in VTK

```
--- a/Utilities/octree/octree/octree_node.txx
+++ b/Utilities/octree/octree/octree_node.txx
@@ -210,7 +210,7 @@ const octree_node<T_, d_, A_>& octree_node<T_, d_, A_>::operator[](int child) co
   {
     throw std::domain_error("Attempt to access children of an octree leaf node.");
   }
-  return this->_M_chilren[child];
+  return this->m_children[child];
 }
```
