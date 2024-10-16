    -p  Tetrahedralizes a piecewise linear complex (PLC). 输入PLC
    -Y  Preserves the input surface mesh (does not modify it).  不在边界上加点
    -r  Reconstructs a previously generated mesh. 重构优化mesh,   tetgen -rq1.2a1e-5V  node,ele,face,edge ,    支持.mesh
    -q  Refines mesh (to improve mesh quality). 提高质量，外接圆半径除以最短边，缺省参数为 2, 1.2, 与角度等价, (1-2) -pq1.2
    -R  Mesh coarsening (to reduce the mesh elements).   # 降低单元数量，-rRV  demo.mesh
    -A  Assigns attributes to tetrahedra in different regions. 不同区域，添加标记，随机顺序，自己对应
    -a  Applies a maximum tetrahedron volume constraint.  全局加密，四面体体积 -a0.001   -a1e-4
    -m  Applies a mesh sizing function.  让.mtr文件生效


    -i  Inserts a list of additional points.
    -O  Specifies the level of mesh optimization.  不做smooth,不移动里面的点
    -s  内部options, 不做smooth,不移动里面的点， tetgen -s0

    -S  Specifies maximum number of added points. 
    -T  Sets a tolerance for coplanar test (default 1e-8). 用于检测共面
    -X  Suppresses use of exact arithmetic.
    -M  No merge of coplanar facets or very close vertices.
    -w  Generates weighted Delaunay (regular) triangulation.
    -c  Retains the convex hull of the PLC.
    -d  Detects self-intersections of facets of the PLC.

    -z  Numbers all output items starting from zero.
    -f  Outputs all faces to .face file.
    -e  Outputs all edges to .edge file.
    -n  Outputs tetrahedra neighbors to .neigh file.
    -v  Outputs Voronoi diagram to files.
    -g  Outputs mesh to .mesh file for viewing by Medit.  输出.mesh
    -k  Outputs mesh to .vtk file for viewing by Paraview.
    -J  No jettison of unused vertices from output .node file.
    -B  Suppresses output of boundary information.
    -N  Suppresses output of .node file. 不要输出node
    -E  Suppresses output of .ele file.  不要输出ele
    -F  Suppresses output of .face and .edge file. 不要输出 face 和 edge
    -I  Suppresses mesh iteration numbers.
    -C  Checks the consistency of the final mesh.
    -Q  Quiet:  No terminal output except errors.
    -V  Verbose:  Detailed information, more terminal output. 输出详细信息，最长边比最小高度， 关键看aspect ratio
    -h  Help:  A brief instruction for using TetGen.
 


增加.ele , .face   .node,  的标记是可以额外增加tag

.mtr文件，人工定义mesh size， 在点上，和点对齐， 楞上现行插值， 绝对的size，   -m 让.mtr生效


情况，
给一个额外的mesh size function 背景网格，比原网格大 a.b.mesh 就是background mesh
