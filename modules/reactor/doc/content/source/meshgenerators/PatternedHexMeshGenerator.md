# PatternedHexMeshGenerator

!syntax description /Mesh/PatternedHexMeshGenerator

## Overview

The `PatternedHexMeshGenerator` object generates a 2D mesh by stitching hexagonal meshes generated by [`PolygonConcentricCircleMeshGenerator`](/PolygonConcentricCircleMeshGenerator.md), [`HexagonConcentricCircleAdaptiveBoundaryMeshGenerator`](/HexagonConcentricCircleAdaptiveBoundaryMeshGenerator.md), and itself, based on a user-defined hexagonal grid pattern. The generated 2D mesh can optionally contain an extra background region and/or external duct regions to form a hexagonal external shape rather than a jagged boundary edge.

## Geometry Information

In order to generate the hexagonal patterned mesh, three fundamental parameters are needed:

- [!param](/Mesh/PatternedHexMeshGenerator/inputs): a vector of mesh generator names that will be used as elements to form the patterned mesh;
- [!param](/Mesh/PatternedHexMeshGenerator/pattern): a two-dimensional vector that represents the pattern of mesh to be generated. The elements must be integers from 0 to N-1, where N is the length of vector parameter [!param](/Mesh/PatternedHexMeshGenerator/inputs).
- [!param](/Mesh/PatternedHexMeshGenerator/pattern_boundary): a MooseEnum parameter that can be either `hexagon` or `none`. When `none` is selected, this object only stitches [!param](/Mesh/PatternedHexMeshGenerator/inputs) meshes into a patterned mesh without extraneous geometry, resulting in an outer boundary with a zig-zag edge. When `hexagon` is selected, a background region is added so that the generated mesh has a hexagonal shape instead of a zig-zag one. Concentric duct regions may also be optionally added to the hexagon periphery when this parameter is `hexagon`. The user can refer to [Figure 1](#pattern_hex) for more details.

!media reactor/meshgenerators/pattern_hex.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=pattern_hex
      caption=A schematic showing the difference between `none` and `hexagon` `pattern_boundary`.

When [!param](/Mesh/PatternedHexMeshGenerator/pattern_boundary) is set as `hexagon`, the user can also provide [!param](/Mesh/PatternedHexMeshGenerator/duct_sizes) in either `apothem` or `radius` style to add external duct regions to the generated hexagon mesh. Also, [!param](/Mesh/PatternedHexMeshGenerator/hexagon_size) must also be provided to define the external hexagon boundary size, which can be either `apothem` or `radius` of the hexagon, as determined by [!param](/Mesh/PatternedHexMeshGenerator/hexagon_size_style). In some cases, if [!param](/Mesh/PatternedHexMeshGenerator/hexagon_size) is small enough, the hexagon external boundary may cut off some of the stitched hexagonal meshes. As long as the concentric circular (`rings`) regions are not cut off, the rest of the mesh is deformed to accommodate the limited external boundaries. Users can also set [!param](/Mesh/PatternedHexMeshGenerator/deform_non_circular_region) as `false` to prevent the non-circular regions from being deformed.

## Control Drum Related MeshMetaData

One of the applications of this object is to generate meshes for prismatic reactor cores. In that case, by setting [!param](/Mesh/PatternedHexMeshGenerator/generate_core_metadata) as true, control drum meshes can also be used as part of [!param](/Mesh/PatternedHexMeshGenerator/inputs) to construct the core mesh. To facilitate the use of control drum rotation simulation objects, a series of `MeshMetaData` can be generated, including:

- `control_drum_positions`: a vector of control drum center positions. This `MeshMetaData` can also be outputted as an ASCII file by setting [!param](/Mesh/PatternedHexMeshGenerator/generate_control_drum_positions_file) as true and providing [!param](/Mesh/PatternedHexMeshGenerator/position_file);
- `control_drum_angles`: a vector of the azimuthal angles of the control drum center positions to the center of the core.
- `control_drums_azimuthal_meta`: a two-dimensional vector containing the sorted azimuthal angles of nodes to the corresponding control drum center for all the control drums.

In addition, [!param](/Mesh/PatternedHexMeshGenerator/assign_control_drum_id) can be set as true so that the control drum [!param](/Mesh/PatternedHexMeshGenerator/inputs) meshes can be indexed using an element extra integer called `control_drum_id`. As illustrated in [Figure 2](#cd_id), the `control_drum_id` is indexed based on the azimuthal angles of the control drums.

!media reactor/functions/azi_cd_id.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=cd_id
      caption=A schematic drawing the indexing rule of `control_drum_id` in the `PatternedHexMeshGenerator` object.

These `MeshMetaData` as well as `control_drum_id` can be used by other MOOSE objects such as [`MultiControlDrumFunction`](/MultiControlDrumFunction.md) to simulate control drums rotation during power transients.

## Interface Boundaries

The user can also decide whether the interface boundaries are generated or not in the peripheral region.

!include /PolygonConcentricCircleMeshGenerator.md start=There are two types end=The user can set

The user can set [!param](/Mesh/PatternedHexMeshGenerator/create_inward_interface_boundaries) and [!param](/Mesh/PatternedHexMeshGenerator/create_outward_interface_boundaries) to control which interface boundaries will be created. If generated, the outward interface boundaries will be assigned ids using sequential odd numbers (i.e., 1, 3, 5, 7, ...) shifted by `INTRINSIC_SIDESET_ID::SLICE_ALT`=30500 from center to periphery, while the inward interface boundaries will be assigned ids using sequential even numbers (i.e., 0, 2, 4, 6, ...) shifted by `INTRINSIC_SIDESET_ID::SLICE_ALT` similarly. 

## Reporting ID Information

This object can generate a hexagonal lattice mesh with `reporting ID` assignments, and can be used successively on its own output mesh to add IDs on the pin and assembly levels, for example.
The reporting ID option can be turned on by defining the name of the reporting ID variable is provided through [!param](/Mesh/PatternedHexMeshGenerator/id_name).

A user can select an ID assignment scheme using [!param](/Mesh/PatternedHexMeshGenerator/assign_type), and the following schemes are currently available:

- `cell` (default):  Assign unique IDs for each component/tile in the lattice in sequential order.

- `pattern`:  Assign IDs based on the ID of the input tiles.

- `manual`: Assign IDs based on user-defined mapping defined in [!param](/Mesh/PatternedHexMeshGenerator/id_pattern).

The default numbering scheme starts at 0 in the upper left hand corner of the hexagon grid (not including duct region) and increments by 1 as the grid is traversed left to right, top to bottom.
In presence of duct regions, separate reporting IDs are automatically generated for the elements in duct regions.
For the `pattern` scheme, all tiles in the pattern with the same input will bear the same reporting ID.
The duct regions will be assigned reporting IDs starting from the next integer higher than the highest one used inside of the ducts.

The name of the reporting ID variable is provided through [!param](/Mesh/PatternedHexMeshGenerator/id_name) depending on the hierarchical level of component.
The ID values themselves are stored as extra element integers on the mesh.
For example, the reporting IDs for individual pins (`pin_id`) can be assigned when assemblies are built because the IDs for pin level are uniquely determined from the pin arrangement within each assembly type.
Similarly, the assembly reporting IDs (`assembly_id`) are assigned in the core construction process.

The multiple `reporting IDs` can be assigned by defining the multiple names of the reporting ID variable, which are provided through the[!param](/Mesh/PatternedHexMeshGenerator/id_name).
The corresponding assignment scheme should be provided in [!param](/Mesh/PatternedHexMeshGenerator/assign_type) for each reporting ID names, accordingly.
In the case that multiple `manual` `assign_type`s are used, the same number of manual ID patterns should be provided in [!param](/Mesh/PatternedHexMeshGenerator/id_pattern).
Each manual pattern in [!param](/Mesh/PatternedHexMeshGenerator/id_pattern) should be separated by using the delimiter '|'.
These defined ID patterns are sequentially assigned to the reporting IDs having `manual` assignment scheme.
The below is an example of using multiple `reporting ID` assignment.
Here, `manual_1_id` uses the first pattern in defined in [!param](/Mesh/PatternedHexMeshGenerator/id_pattern), and `manual_2_id` uses the second one.

!listing!
id_name = 'manual_1_id cell_id manual_2_id'
assign_type 'manual cell manual'
id_pattern = '1 1;
             2 2 2;
              3 3|
              1 2;
             1 2 3
              2 3;
!listing-end!

Certain regions can be excluded from being labeled with an ID, for example dummy regions that will later be deleted.
This can be accommodated by listing mesh objects in the [!param](/Mesh/PatternedHexMeshGenerator/exclude_id) input parameter.
IDs will not be assigned to these mesh objects.
Usage of this parameter is helpful to retain sequential numbering when dummy region are later deleted, or to only label areas of interest.


## Example Syntax

!listing modules/reactor/test/tests/meshgenerators/patterned_hex_mesh_generator/patterned_pattern.i block=Mesh/pattern_1

!syntax parameters /Mesh/PatternedHexMeshGenerator

!syntax inputs /Mesh/PatternedHexMeshGenerator

!syntax children /Mesh/PatternedHexMeshGenerator
