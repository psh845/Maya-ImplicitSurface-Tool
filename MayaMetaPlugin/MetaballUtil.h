#include <maya/MStatus.h>
#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MFloatArray.h>
#include <maya/MTypes.h>
#include <maya/MGlobal.h>
#include <maya/MVector.h>
#include <maya/MFnMesh.h>
#include <maya/MString.h>
#include <maya/MItSelectionList.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MDagPath.h>
#include <maya/MFnTransform.h>
#include <maya/MSelectionList.h>
#include <maya/MItMeshEdge.h>
#include <maya/MItDag.h>
#include <maya/MPlug.h>
#include <maya/MFnDependencyNode.h>
#include <math.h>



#include<time.h> // time() 함수 사용위함

class GridCell{
public:
	MPointArray p;
	MFloatArray val;
	MColorArray c;
	MColor cc;
	int add_C;

	GridCell()
	{
		p.clear();
		val.clear();
		c.clear();
		cc=(0.0,0.0,0.0);
		add_C=0;
	};

	~GridCell()
	{

	};

	void set(MPoint point, float density )
	{
		
		p.append(point);
		val.append(density);
	};
	
	void setColor(MColor C)
	{
		
		c.append(C);
		add_C++;
		cc.operator+=(C);
	};

	void divColor()
	{
		cc.operator/=(add_C);
	}

	
};

MStatus Cube_Array(
	float &Xarray,
	float &Yarray,
	float &Zzrray,
	float &width,
	MFloatArray &XArrayList,
	MFloatArray &YArrayList,
	MFloatArray &ZArrayList
	);

MStatus Cube_box(
	MPointArray &ArrayVertices,
	MFloatArray &XArrayList,
	MFloatArray &YArrayList,
	MFloatArray &ZArrayList,
	MIntArray &isoCheckList
	);

MStatus sBox();

float Metaball(
	MPoint &boxVertice, 
	float &width,
	MPointArray &centerSphere,
	MFloatArray &radius,
	MFloatArray &weight,
	MColorArray &color,
	float &isolevel,
	MColor &c
	);

float Metaball_line(
	MPoint &boxVertice, 
	float &width,
	MPointArray &line_point,
	MFloatArray &line_radius,
	MFloatArray &line_radius2,
	MFloatArray &line_weight
	);

float Metaball_curve4p(
	MPoint &boxVertice,
	float &width,
	MPointArray &curve_point,
	MFloatArray &curve_radius,
	MFloatArray &curve_radius2,
	MFloatArray &curve_weight
	);

float Metaball_curve(
	MPoint &boxVertice,
	float &width,
	MPointArray &curve_point,
	MFloatArray &curve_radius,
	MFloatArray &curve_radius2,
	MFloatArray &curve_weight,
    MIntArray &curve_npoint
	);

float Metaball_curveNurbs(
	MSelectionList &selection,
	MPoint &boxVertice,
	float &width,
	MFloatArray &curve_radius,
	MFloatArray &curve_radius2,
	MFloatArray &curve_weight,
	MColorArray &curve_color,
	float &isolevel,
	MColor &c
	);

MStatus MarchingCube(
	GridCell &grid,
	float &isolevel,
	MPointArray &triangles,
	MIntArray &polyCount,
	MIntArray &polyconnet,
	int &polyconnetCount,
	MColorArray &trianglesColor
	);

MPoint VertexInterp(float &isolevel,MPoint &p1,MPoint &p2, float &valp1, float  &valp2);
MColor ColorInterp(float isolevel, MColor c1, MColor c2, float valp1, float valp2);

void ColorMeta(MColor &c,float w, MColor Curve_color,float isolevel);

bool gridCheck(
	MIntArray &isoCheckList,
	int val);

MStatus getData_curve(
	MSelectionList &selection,
	MFloatArray &curve_radius,
	MFloatArray &curve_radius2,
	MFloatArray &curve_weight,
	MColorArray &curve_color
	);


MStatus getData(
	MSelectionList &selection,
	MPointArray &centerSphere,
	MFloatArray &radius,
	MFloatArray &weight,
	MColorArray &color
	);