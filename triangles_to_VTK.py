#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from evtk.vtk import VtkFile, VtkUnstructuredGrid, VtkTriangle, VtkVertex
import numpy as np

from copy import deepcopy


"""
UnstructuredGrid
================
Each 'UnstructuredGrid' piece specifies a set of points and cells independently
from the other pieces. The points are described explicitly by the 'Points'
element. The cells are described explicitly by the 'Cells' element.

    <VTKFile type=”UnstructuredGrid” ...>
        <UnstructuredGrid>
            <Piece NumberOfPoints=”#” NumberOfCells=”#”>
                <PointData>...</PointData>
                <CellData>...</CellData>
                <Points>...</Points>
                <Cells>...</Cells>
            </Piece>
        </UnstructuredGrid>
    </VTKFile>

Every dataset describes the data associated with its points and cells with
'PointData' and 'CellData' XML elements as follows:

    <PointData Scalars=”Temperature” Vectors=”Velocity”>
        <DataArray Name=”Velocity” .../>
        <DataArray Name=”Temperature” .../>
        <DataArray Name=”Pressure” .../>
    </PointData>

VTK allows an arbitrary number of data arrays to be associated with the points
and cells of a dataset. Each data array is described by a 'DataArray' element
which, among other things, gives each array a name. The following attributes of
'PointData' and 'CellData' are used to specify the active arrays by name:
    'Scalars' — The name of the active scalars array, if any.
    'Vectors' — The name of the active vectors array, if any.
    'Normals' — The name of the active normals array, if any.
    'Tensors' — The name of the active tensors array, if any.
    'TCoords' — The name of the active texture coordinates array, if any.

Points
======
The 'Points' element explicitly defines coordinates for each point
individually. It contains one 'DataArray' element describing an array with
three components per value, each specifying the coordinates of one point.
    <Points>
        <DataArray NumberOfComponents=”3” .../>
    </Points>

Coordinates
===========
The 'Coordinates' element defines point coordinates for an extent by specifying
the ordinate along each axis for each integer value in the extent’s range. It
contains three 'DataArray' elements describing the ordinates along the
x-y-z axes, respectively.
    <Coordinates>
        <DataArray .../>
        <DataArray .../>
        <DataArray .../>
    </Coordinates>

Verts, Lines, Strips, and Polys
===============================
The 'Verts', 'Lines', 'Strips', and 'Polys' elements define cells explicitly by
specifying point connectivity. Cell types are implicitly known by the type of
element in which they are specified. Each element contains two 'DataArray'
elements. The first array specifies the point connectivity. All the cells’
point lists are concatenated together. The second array specifies the offset
into the connectivity array for the end of each cell.
    <Verts>
        <DataArray type=”Int32” Name=”connectivity” .../>
        <DataArray type=”Int32” Name=”offsets” .../>
    </Verts>

Cells
=====
The 'Cells' element defines cells explicitly by specifying point connectivity
and cell types. It contains three 'DataArray' elements. The first array
specifies the point connectivity. All the cells’ point lists are concatenated
together. The second array specifies the offset into the connectivity array
for the end of each cell. The third array specifies the type of each cell.
(Note: the cell types are defined in Figure 2 and Figure 3.)
    <Cells>
        <DataArray type=”Int32” Name=”connectivity” .../>
        <DataArray type=”Int32” Name=”offsets” .../>
        <DataArray type=”UInt8” Name=”types” .../>
    </Cells>

All of the data and geometry specifications use 'DataArray' elements to
describe their actual content as follows:

DataArray
=========
The 'DataArray' element stores a sequence of values of one type. There may be
one or more components per value.
    <DataArray type=”Float32” Name=”vectors” NumberOfComponents=”3”
               format=”appended” offset=”0”/>
    <DataArray type=”Float32” Name=”scalars” format=”binary”>
    bAAAAAAAAAAAAIA/AAAAQAAAQEAAAIBA... </DataArray>
    <DataArray type=”Int32” Name=”offsets” format=”ascii”>
    10 20 30 ... </DataArray>

The attributes of the 'DataArray' elements are described as follows:
    type — The data type of a single component of the array. This is one of
           Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32,
           Float64.

           Note: the 64-bit integer types are only supported if
           VTK_USE_64BIT_IDS is on (a CMake variable—see “CMake” on page 8) or
           the platform is 64-bit.

    Name — The name of the array. This is usually a brief description of the
           data stored in the array.

    NumberOfComponents — The number of components per value in the array.

    format — The means by which the data values themselves are stored in the
             file. This is “ascii”, “binary”, or “appended”.

    offset — If the format attribute is “appended”, this specifies the offset
             from the beginning of the appended data section to the beginning
             of this array’s data.

The format attribute chooses among the three ways in which data values can be
stored:
    format=”ascii” — The data are listed in ASCII directly inside the DataArray
                     element. Whitespace is used for separation.
    format=”binary” — The data are encoded in base64 and listed contiguously
                      inside the DataArray element.
                      Data may also be compressed before encoding in base64.
                      The byte-order of the data matches that specified by
                      the byte_order attribute of the 'VTKFile' element.
    format=”appended” — The data are stored in the appended data section.
                        Since many 'DataArray' elements may store their data in
                        this section, the 'offset' attribute is used to specify
                        where each DataArray’s data begins.

This format ("appended") is the default used by VTK’s writers.

The appended data section is stored in an 'AppendedData' element that is nested
inside 'VTKFile' after the dataset element:
<VTKFile ...>
...
    <AppendedData encoding=”base64”>
        _QMwEAAAAAAAAA...
    </AppendedData>
</VTKFile>

The appended data section begins with the first character after the underscore
inside the AppendedData element. The underscore is not part of the data, but is
always present. Data in this section is always in binary form, but can be
compressed and/or base64 encoded. The byte-order of the data matches that
specified by the byte_order attribute of the 'VTKFile' element. Each
DataArray’s data are stored contiguously and appended immediately after the
previous DataArray’s data without a separator. The DataArray’s offset attribute
indicates the file position offset from the first character after the
underscore to the beginning its data.
"""


# =================================
#       Helper functions
# =================================
def _addDataToFile(vtkFile, cellData, pointData):
    # Point data
    if pointData:
        keys = sorted(list(pointData.keys()))
        vtkFile.openData("Point", scalars = keys[0])
        for key in keys:
            data = pointData[key]
            vtkFile.addData(key, data)
        vtkFile.closeData("Point")

    # Cell data
    if cellData:
        keys = sorted(list(cellData.keys()))
        vtkFile.openData("Cell", scalars = keys[0])
        for key in keys:
            data = cellData[key]
            vtkFile.addData(key, data)
        vtkFile.closeData("Cell")

def _appendDataToFile(vtkFile, cellData, pointData):
    # Append data to binary section
    if pointData != None:
        keys = sorted(list(pointData.keys()))
        for key in keys:
            data = pointData[key]
            vtkFile.appendData(data)

    if cellData != None:
        keys = sorted(list(cellData.keys()))
        for key in keys:
            data = cellData[key]
            vtkFile.appendData(data)
            
   

def triangle_faces_to_VTK(filename, x, y, z, faces, point_data=None, cell_data=None):
    vertices = (x, y, z)

    w = VtkFile(filename, VtkUnstructuredGrid)
    w.openGrid()
    w.openPiece(npoints=len(x), ncells=len(faces))
    w.openElement("Points")
    w.addData("Points", vertices)
    w.closeElement("Points")

    # Create some temporary arrays to write grid topology.
    ncells = len(faces)
    # Index of last node in each cell.
    offsets = np.arange(start=3, stop=3*(ncells + 1), step=3, dtype='uint32')
    # Connectivity as unrolled array.
    connectivity = faces.reshape(ncells*3).astype('int32')
    cell_types = np.ones(ncells, dtype='uint8')*VtkTriangle.tid

    w.openElement("Cells")
    w.addData("connectivity", connectivity)
    w.addData("offsets", offsets)
    w.addData("types", cell_types)
    w.closeElement("Cells")

    _addDataToFile(w, cellData=cell_data, pointData=point_data)

    w.closePiece()
    w.closeGrid()

    w.appendData(vertices)
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)

    _appendDataToFile(w, cellData=cell_data, pointData=point_data)

    w.save()
    return w.getFileName()


def pointsToVTK(path, x, y, z, data):
    """
        Export points and associated data as an unstructured grid.

        PARAMETERS:
            path: name of the file without extension where data to be saved.
            x, y, z: 1D arrays with coordinates of the points.
            data: dictionary {'varname': data_array} of point-data to export.

        RETURNS:
            Full path to saved file.
    """
    assert (x.size == y.size == z.size)
    npoints = x.size

    # Create some temporary arrays to write grid topology ...
    # ... index of last node in each cell ...
    offsets = np.arange(start=1, stop=npoints + 1, dtype='int32')
    # ... unwinding our triangles into connectivity list.
    connectivity = np.arange(npoints, dtype='int32')
    cell_types = np.empty(npoints, dtype='uint8')
    cell_types[:] = VtkVertex.tid

    w = VtkFile(path, VtkUnstructuredGrid)
    w.openGrid()
    w.openPiece(ncells=npoints, npoints=npoints)

    w.openElement("Points")
    w.addData("points", (x, y, z))
    w.closeElement("Points")
    w.openElement("Cells")
    w.addData("connectivity", connectivity)
    w.addData("offsets", offsets)
    w.addData("types", cell_types)
    w.closeElement("Cells")

    _addDataToFile(w, cellData=None, pointData=data)

    w.closePiece()
    w.closeGrid()
    w.appendData((x, y, z))
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)

    _appendDataToFile(w, cellData=None, pointData=data)

    w.save()
    return w.getFileName()


class Level:
    """Approximation of the spherical surface made from triangles.

    It support refinement of the surface to increase detalization uniformly by
    decreasing face size 4 times on each step. Enumeration of vertices within
    face is shown below together with the numbers of new vertices introduced
    during refinement:

                2                                           2
                *                                           *
               / \                       .                 / \
              /   \                      |\               /   \
             /     \           +---------+ \           5 /     \ 4
            /       \          |            *           *-------*
           /         \         +---------+ /           / \     / \
          /           \                  |/           /   \   /   \
         /             \                 "           /     \ /     \
        *---------------*                           *-------*-------*
       0                 1                         0        3        1

    """
    def __init__(self, faces, vertices):
        self.faces = faces
        self.vertices = vertices

    def refine(self):
        """Returns level with better sampling."""
        cache = {}
        vert = deepcopy(self.vertices)

        def _(i, j):
            '''Makes cached middle vertex and returns its index.'''
            ij = tuple(sorted([i, j]))
            if ij not in cache:
                pos = len(vert)
                v = 0.5*(vert[i] + vert[j])
                vert.append(v/np.sqrt(np.inner(v, v)))
                cache[ij] = pos
            return cache[ij]

        f_pos = 0
        face = np.zeros([4*len(self.faces), 3], int)
        for f in self.faces:
            i, j, k = f
            face[f_pos+0, :] = [  i,     _(i, j), _(i, k)]
            face[f_pos+1, :] = [_(i, j),   j,     _(j, k)]
            face[f_pos+2, :] = [_(i, j), _(j, k), _(i, k)]
            face[f_pos+3, :] = [  k,     _(i, k), _(j, k)]
            f_pos += 4
        return Level(face, vert)


class Polyhedron:
    """http://student.ulb.ac.be/~claugero/sphere/index.html"""
    def __init__(self, level=3, base="tetrahedron"):
        """Creates set of normals."""
        assert level in range(1, 8), "bad level: %s" % level

        if base == "tetrahedron":
            s = 1.0/np.sqrt(3.0);
            verts = [[ s,  s,  s], [-s, -s,  s],
                     [-s,  s, -s], [ s, -s, -s]]
            faces = [[0, 2, 1], [0, 1, 3], [2, 3, 1], [3, 2, 0]]
        elif base == "octahedron":
            verts = [[ 0.0,  0.0, -1.0], [ 1.0,  0.0,  0.0],
                     [ 0.0, -1.0,  0.0], [-1.0,  0.0,  0.0],
                     [ 0.0,  1.0,  0.0], [ 0.0,  0.0,  1.0] ]
            faces = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1],
                     [5, 2, 1], [5, 3, 2], [5, 4, 3], [5, 1, 4]]
        elif base == "icosahedron":
            t = (1 + np.sqrt(5))/2
            tau = t/np.sqrt(1 + t*t)
            one = 1/np.sqrt(1 + t*t)
            verts = [[ tau,  one,  0.0], [-tau,  one,  0.0],
                     [-tau, -one,  0.0], [ tau, -one,  0.0],
                     [ one,  0.0,  tau], [ one,  0.0, -tau],
                     [-one,  0.0, -tau], [-one,  0.0,  tau],
                     [ 0.0,  tau,  one], [ 0.0, -tau,  one],
                     [ 0.0, -tau, -one], [ 0.0,  tau, -one]]

            faces = [[4, 8, 7], [4, 7, 9], [5, 6, 11], [5, 10, 6],
                     [0, 4, 3], [0, 3, 5], [2, 7, 1], [2, 1, 6],
                     [8, 0, 11], [8, 11, 1], [9, 10, 3], [9, 2, 10],
                     [8, 4, 0], [11, 0, 5], [4, 9, 3], [5, 3, 10],
                     [7, 8, 1], [6, 1, 11], [7, 2, 9], [6, 10, 2]]
        else:
            raise ValueError("unknown base figure: " + str(base))

        self.level = level
        self.base = base
        self.tree = [Level(faces=faces, vertices=[np.array(v) for v in verts])]
        for l in range(1, level + 1):
            self.tree.append(self.tree[-1].refine())

    def get_approximation(self, n):
        return self.tree[n].faces, self.tree[n].vertices


# import ipdb; ipdb.set_trace()

# arr = np.array
# sphere = Polyhedron(level=3, base="icosahedron")
# faces, vertices = sphere.get_approximation(3)

# x = arr([v[0] for v in vertices])
# y = arr([v[1] for v in vertices])
# z = arr([v[2] for v in vertices])
# cell_scalar = arr([x[i] + x[j] + x[k] for i, j, k in faces])
# triangle_faces_to_VTK("demo_sphere",
#                       x=x, y=y, z=z,
#                       faces=arr(faces),
#                       point_data={'s': x*x + y*y,
#                                   'n': (x, y, z),
#                                   'v': (-y, x, -np.sin(np.arccos(-1.0)*z))},
#                       cell_data={'rho': cell_scalar})

# pointsToVTK("demo_stock",
#             x, y, z,
#             {"rho": x*x + y*y})
