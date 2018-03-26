//
// Copyright (C) 
// 
// File: MayaMetaPluginCmd.cpp
//
// MEL Command: MayaMetaPlugin
//
// Author: Maya Plug-in Wizard 2.0
//

// Includes everything needed to register a simple MEL command with Maya.
// 
#include <maya/MSimple.h>
#include <maya/MArgList.h>
#include <maya/MPxCommand.h>
#include <maya/MFnPlugin.h>
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MDGModifier.h>
#include <maya/MProgressWindow.h>
//#include <maya/MSelectionList.h>
//#include <maya/MGlobal.h>
//#include <maya/MFnMesh.h>
//#include <maya/MFnNurbsCurve.h>
//#include <maya/MFnNurbsSurface.h>
//#include <maya/MDagPath.h>
//#include <maya/MPointArray.h>
//#include <maya/MItSelectionList.h>
//#include <maya/MItMeshEdge.h>
#include <time.h>
#include <iostream>
#include "MetaballUtil.h"

// Use helper macro to register a command with Maya.  It creates and
// registers a command that does not support undo or redo.  The 
// created class derives off of MPxCommand.
//
//DeclareSimpleCommand( MayaMetaPlugin, "", "2014");

class MayaMetaPlugin : public MPxCommand 
{
public:
	virtual MStatus	doIt ( const MArgList& );
	virtual MStatus undoIt();
 	virtual MStatus redoIt();
	virtual bool isUndoable() const { return true; } 

	static void *creator() { return new MayaMetaPlugin; }
	static MSyntax newSyntax();
private:
	MDGModifier dgMod;
};


const char *sizeXFlag = "-sx", *sizeXLongFlag = "-sizeX";
const char *sizeYFlag = "-sy", *sizeYLongFlag = "-sizeY";
const char *sizeZFlag = "-sz", *sizeZLongFlag = "-sizeZ";
const char *resolutionFlag = "-rs", *resolutionLongFlag = "-resolution";
const char *thresholdFlag = "-th", *thresholdLongFlag = "-threshold";


MSyntax MayaMetaPlugin::newSyntax()
{
    MSyntax syntax;

	syntax.addFlag( sizeXFlag, sizeXLongFlag, MSyntax::kDouble );
    syntax.addFlag( sizeYFlag, sizeYLongFlag, MSyntax::kDouble );
	syntax.addFlag( sizeZFlag, sizeZLongFlag, MSyntax::kDouble );
	syntax.addFlag( resolutionFlag, resolutionLongFlag, MSyntax::kDouble );
	syntax.addFlag( thresholdFlag, thresholdLongFlag, MSyntax::kDouble );

	return syntax;
}

const char *helpText =
"\nFor further details consult the help documentation."
"\nFor quick help use: help MayaMetaPlugin";

MStatus MayaMetaPlugin::doIt( const MArgList &args )
//
//	Description:
//		implements the MEL MayaMetaPlugin command.
//
//	Arguments:
//		args - the argument list that was passes to the command from MEL
//
//	Return Value:
//		MS::kSuccess - command succeeded
//		MS::kFailure - command failed (returning this value will cause the 
//                     MEL script that is being run to terminate unless the
//                     error is caught using a "catch" statement.
//
{
	clock_t begin,end;
	begin = clock();

	float Xarry = 20.0;
	float Yarry = 20.0;
	float Zarry = 20.0;
	float width = 1.0;
	float threshold = 1.0;
	int polyconnetCount =0;

	double _Xarry=0.0;
	double _Yarry=0.0;
	double _Zarry=0.0;
	double _width=0.0;
	double _threshold=0.0;


	MArgDatabase argData( syntax(), args );

	
	if( argData.isFlagSet( sizeXFlag ) )
		argData.getFlagArgument( sizeXFlag, 0, _Xarry );


	if( argData.isFlagSet( sizeYFlag ) )
		argData.getFlagArgument( sizeYFlag, 0, _Yarry );

	if( argData.isFlagSet( sizeZFlag ) )
		argData.getFlagArgument( sizeZFlag, 0, _Zarry );

	if( argData.isFlagSet( resolutionFlag ) )
		argData.getFlagArgument( resolutionFlag, 0, _width );

	if( argData.isFlagSet( thresholdFlag ) )
		argData.getFlagArgument( thresholdFlag, 0, _threshold );

	if(_Xarry>0)
		Xarry = _Xarry;
	if(_Yarry>0)
		Yarry = _Yarry;
	if(_Zarry>0)
		Zarry = _Zarry;
	if(_width>0)
		width = _width;
	if(_threshold>0)
		threshold = _threshold;


	MStatus stat = MS::kSuccess;

	MFloatArray XArrayList;
	MFloatArray YArrayList;
	MFloatArray ZArrayList;

	MFloatArray isoValueList;
	MIntArray isoCheckList;
	MColorArray isoColorList;
	

	XArrayList.clear();
	YArrayList.clear();
	ZArrayList.clear();

	isoValueList.clear();
	isoColorList.clear();
	isoCheckList.clear();


	MSelectionList selection;
	selection.clear();
	MGlobal::getActiveSelectionList(selection);
	

	//MPointArray boxVertices;
	MPointArray ArrayVertices;

	MPointArray centerSphere;
	MFloatArray radius;
	MFloatArray weight;
	MColorArray color;

	MPointArray line_point;
	MFloatArray line_radius;
	MFloatArray line_radius2;
	MFloatArray line_weight;

	MPointArray Curve_point;
	MFloatArray Curve_radius;
	MFloatArray Curve_radius2;
	MFloatArray Curve_weight;
	MIntArray Curve_npoint;
	MColorArray Curve_color;

	MPointArray MC;
	MFloatArray IsoValue;

	MPointArray triangles;
	MIntArray polyCount;
	MIntArray polyconnet;
	MColorArray trianglesColor;

	//boxVertices.clear();
	ArrayVertices.clear();
	IsoValue.clear();
	centerSphere.clear();
	radius.clear();
	weight.clear();
	color.clear();

	line_point.clear();
	line_radius.clear();
	line_radius2.clear();
	line_weight.clear();
	
	Curve_point.clear();
	Curve_radius.clear();
	Curve_radius2.clear();
	Curve_weight.clear();
	Curve_npoint.clear();
	Curve_color.clear();

	MC.clear();

	IsoValue.setLength(8);
	MC.setLength(8);

	triangles.clear();
	polyCount.clear();
	polyconnet.clear();
	trianglesColor.clear();
	
	/*
	centerSphere.append(MPoint(5.0,3.0,0.0));
	//centerSphere.append(MPoint(1.2,2.5,0.0));
	radius.append(1.0);
	//radius.append(2.0);
	weight.append(2.0);
	//weight.append(2.0);

	line_point.append(MPoint(1.2,-1.0,1.0));
	line_point.append(MPoint(1.2,-4.0,1.0));
	line_radius.append(1.0);
	line_weight.append(2.0);

	line_point.append(MPoint(1.2,-1.0,-1.0));
	line_point.append(MPoint(1.2,-4.0,-1.0));
	line_radius.append(1.0);
	line_weight.append(2.0);

	line_point.append(MPoint(-1.0,-1.0,1.0));
	line_point.append(MPoint(-1.0,-4.0,1.0));
	line_radius.append(1.0);
	line_weight.append(2.0);

	line_point.append(MPoint(-1.0,-1.0,-1.0));
	line_point.append(MPoint(-1.0,-4.0,-1.0));
	line_radius.append(1.0);
	line_weight.append(2.0);

	line_point.append(MPoint(1.5,0.0,0.0));
	line_point.append(MPoint(-2.0,0.0,0.0));
	line_radius.append(3.0);
	line_weight.append(2.0);

	Curve_point.append(MPoint(1.5,0.0,0.0));
	Curve_point.append(MPoint(2.6,1.0,0.0));
	Curve_point.append(MPoint(3.2,2.0,0.0));
	Curve_point.append(MPoint(4.0,3.0,0.0));
	Curve_radius.append(1.5);
	Curve_weight.append(2.0);

	Curve_point.append(MPoint(-2,0.0,0.0));
	Curve_point.append(MPoint(-2.5,-0.5,0.0));
	Curve_point.append(MPoint(-3.0,-2.0,0.0));
	Curve_point.append(MPoint(-5.0,-3.0,0.0));
	Curve_radius.append(1.2);
	Curve_weight.append(2.0);*/

	//MItSelectionList iter(selection);

	MDagPath dagPath;
	MObject component;
	MFnNurbsCurve curvFn;
	MFnNurbsSurface surfaceFn;
	MFnTransform transformFn;
	MFnDagNode nodeFn;
	MString name;
	
	stat=getData(selection,centerSphere,radius,weight,color);
	stat=getData_curve(selection,Curve_radius,Curve_radius2,Curve_weight,Curve_color);

	//iter.setFilter(MFn::kNurbsCurve);
	int tt=0;

	//for (;!iter.isDone();iter.next())
	//{
	//	
	//	MString txt;
	//	txt =  MString(" ");
	//	iter.getDagPath(dagPath);
	//	curvFn.setObject(dagPath);

	//	double tStart, tEnd;
	//	

	//	//print(curvFn.getCVs());
	//	int aa =curvFn.numCVs(&stat);
	//	Curve_npoint.append(aa);

	//	curvFn.getKnotDomain(tStart, tEnd);
	//	//stat =curvFn.getCVs(Curve_point,MSpace::kWorld);

	//	Curve_radius.append(2);
	//	Curve_weight.append(2);

	//	
	//	for(int i=0; i<aa; i++)
	//	{
	//		MPoint pt;
	//		curvFn.getCV(i,pt,MSpace::kWorld);
	//		Curve_point.append(pt);
	//		txt += MString(" ") + Curve_point[tt+i].x+Curve_point[tt+i].y+Curve_point[tt+i].z+"\t";

	//	}

	//	
	//	MGlobal::displayInfo( txt );
	//	tt=aa;
	//
	//}

	
	stat = Cube_Array(Xarry,Yarry,Zarry,width,XArrayList,YArrayList,ZArrayList);


	if(!stat)
		stat.perror("Unble to create mesh");
	
	

	for(int x=0; x<XArrayList.length()-1; x++)
	{
		for(int y=0; y<YArrayList.length()-1; y++)
		{
			for(int z=0; z<ZArrayList.length()-1; z++)
			{


				MPoint p;
				
				float sum;
				MColor c;
				
				sum=0;
				c.r=0.0;
				c.g=0.0;
				c.b=0.0;
				p= MPoint(XArrayList[x],YArrayList[y],ZArrayList[z]);
				ArrayVertices.append(p);
				sum = Metaball(p,width,centerSphere,radius,weight,color,threshold,c);
				sum += Metaball_line(p,width,line_point,line_radius,line_radius2,line_weight);
				//sum += Metaball_curve(p,width,Curve_point,Curve_radius,Curve_weight,Curve_npoint);
				sum += Metaball_curveNurbs(selection,p,width,Curve_radius,Curve_radius2, Curve_weight,Curve_color,threshold,c);
				isoValueList.append(sum);
				if(sum>=threshold)
					isoCheckList.append(1);
				else
					isoCheckList.append(0);
				isoColorList.append(c);


				sum=0;
				c.r=0.0;
				c.g=0.0;
				c.b=0.0;
				p= MPoint(XArrayList[x+1],YArrayList[y],ZArrayList[z]);
				ArrayVertices.append(p);
				sum = Metaball(p,width,centerSphere,radius,weight,color,threshold,c);
				sum += Metaball_line(p,width,line_point,line_radius,line_radius2,line_weight);
				//sum += Metaball_curve(p,width,Curve_point,Curve_radius,Curve_weight,Curve_npoint);
				sum += Metaball_curveNurbs(selection,p,width,Curve_radius,Curve_radius2, Curve_weight,Curve_color,threshold,c);
				isoValueList.append(sum);
				if(sum>=threshold)
					isoCheckList.append(1);
				else
					isoCheckList.append(0);
				isoColorList.append(c);


				sum=0;
				c.r=0.0;
				c.g=0.0;
				c.b=0.0;
				p= MPoint(XArrayList[x+1],YArrayList[y],ZArrayList[z+1]);
				ArrayVertices.append(p);
				sum = Metaball(p,width,centerSphere,radius,weight,color,threshold,c);
				sum += Metaball_line(p,width,line_point,line_radius,line_radius2, line_weight);
				//sum += Metaball_curve(p,width,Curve_point,Curve_radius,Curve_weight,Curve_npoint);
				sum += Metaball_curveNurbs(selection,p,width,Curve_radius,Curve_radius2, Curve_weight,Curve_color,threshold,c);
				isoValueList.append(sum);
				if(sum>=threshold)
					isoCheckList.append(1);
				else
					isoCheckList.append(0);
				isoColorList.append(c);
				

				sum=0;
				c.r=0.0;
				c.g=0.0;
				c.b=0.0;
				p= MPoint(XArrayList[x],YArrayList[y],ZArrayList[z+1]);
				ArrayVertices.append(p);
				sum = Metaball(p,width,centerSphere,radius,weight,color,threshold,c);
				sum += Metaball_line(p,width,line_point,line_radius,line_radius2,line_weight);
				//sum += Metaball_curve(p,width,Curve_point,Curve_radius,Curve_weight,Curve_npoint);
				sum += Metaball_curveNurbs(selection,p,width,Curve_radius,Curve_radius2, Curve_weight,Curve_color,threshold,c);
				isoValueList.append(sum);
				if(sum>=threshold)
					isoCheckList.append(1);
				else
					isoCheckList.append(0);
				isoColorList.append(c);


				sum=0;
				c.r=0.0;
				c.g=0.0;
				c.b=0.0;
				p= MPoint(XArrayList[x],YArrayList[y+1],ZArrayList[z]);
				ArrayVertices.append(p);
				sum = Metaball(p,width,centerSphere,radius,weight,color,threshold,c);
				sum += Metaball_line(p,width,line_point,line_radius,line_radius2,line_weight);
				//sum += Metaball_curve(p,width,Curve_point,Curve_radius,Curve_weight,Curve_npoint);
				sum += Metaball_curveNurbs(selection,p,width,Curve_radius,Curve_radius2,Curve_weight,Curve_color,threshold,c);
				isoValueList.append(sum);
				if(sum>=threshold)
					isoCheckList.append(1);
				else
					isoCheckList.append(0);
				isoColorList.append(c);


				sum=0;
				c.r=0.0;
				c.g=0.0;
				c.b=0.0;
				p= MPoint(XArrayList[x+1],YArrayList[y+1],ZArrayList[z]);
				ArrayVertices.append(p);
				sum = Metaball(p,width,centerSphere,radius,weight,color,threshold,c);
				sum += Metaball_line(p,width,line_point,line_radius,line_radius2,line_weight);
				//sum += Metaball_curve(p,width,Curve_point,Curve_radius,Curve_weight,Curve_npoint);
				sum += Metaball_curveNurbs(selection,p,width,Curve_radius,Curve_radius2,Curve_weight,Curve_color,threshold,c);
				isoValueList.append(sum);
				if(sum>=threshold)
					isoCheckList.append(1);
				else
					isoCheckList.append(0);
				isoColorList.append(c);


				sum=0;
				c.r=0.0;
				c.g=0.0;
				c.b=0.0;
				p= MPoint(XArrayList[x+1],YArrayList[y+1],ZArrayList[z+1]);
				ArrayVertices.append(p);
				sum = Metaball(p,width,centerSphere,radius,weight,color,threshold,c);
				sum += Metaball_line(p,width,line_point,line_radius,line_radius2,line_weight);
				//sum += Metaball_curve(p,width,Curve_point,Curve_radius,Curve_weight,Curve_npoint);
				sum += Metaball_curveNurbs(selection,p,width,Curve_radius,Curve_radius2,Curve_weight,Curve_color,threshold,c);
				isoValueList.append(sum);
				if(sum>=threshold)
					isoCheckList.append(1);
				else
					isoCheckList.append(0);
				isoColorList.append(c);


				sum=0;
				c.r=0.0;
				c.g=0.0;
				c.b=0.0;
				p= MPoint(XArrayList[x],YArrayList[y+1],ZArrayList[z+1]);
				ArrayVertices.append(p);
				sum = Metaball(p,width,centerSphere,radius,weight,color,threshold,c);
				sum += Metaball_line(p,width,line_point,line_radius,line_radius2,line_weight);
				//sum += Metaball_curve(p,width,Curve_point,Curve_radius,Curve_weight,Curve_npoint);
				sum += Metaball_curveNurbs(selection,p,width,Curve_radius,Curve_radius2,Curve_weight,Curve_color,threshold,c);
				isoValueList.append(sum);
				if(sum>=threshold)
					isoCheckList.append(1);
				else
					isoCheckList.append(0);
				isoColorList.append(c);

				
			}
		}
	}


	////////////////////////////////////////////////////////////////
	//MString title = "Converting ";
 //    MString sleeping = "waiting: ";
 //       
 //    int amount = 0;
 //    int maxProgress = 100;

 // // Set up and print progress window state

	// MString progressWindowState = MString("Progress Window Info:") +
 //               MString("\nMin: ") + MProgressWindow::progressMin() +
 //               MString("\nMax: ") + MProgressWindow::progressMax() + 
 //               MString("\nTitle: ") + MProgressWindow::title() + 
 //               MString("\nInterruptible: ") + MProgressWindow::isInterruptable();

    /* MGlobal::displayInfo(progressWindowState);

	 CHECK_MSTATUS(MProgressWindow::startProgress());

	 for (int i = amount; i < maxProgress; i++)
        {
                if (i != 0 && MProgressWindow::isCancelled()) 
				{
                        MGlobal::displayInfo("Progress interrupted!");
                        break;
                
				}

	
				MString statusStr = sleeping;
				statusStr += i;
				 CHECK_MSTATUS(MProgressWindow::setProgressStatus(statusStr));
                CHECK_MSTATUS(MProgressWindow::advanceProgress(1));

                MGlobal::displayInfo(MString("Current progress: ") + MProgressWindow::progress());

                MGlobal::executeCommand("pause -sec 1", false, false);
	 }

	CHECK_MSTATUS(MProgressWindow::endProgress());*/

	////////////////////////////////////////////////////////////

	stat = Cube_box(ArrayVertices,XArrayList,YArrayList,ZArrayList,isoCheckList);
	

	if(!stat)
		stat.perror("Unble to create mesh");

	for(int i=0; i <ArrayVertices.length(); i+=8)
	{

		//if(gridCheck(isoCheckList,i))
		//{
			GridCell grid;
			MC.clear();
			IsoValue.clear();


			MC[0]= ArrayVertices[i];
			IsoValue[0] = isoValueList[i];
			grid.set(MC[0],IsoValue[0]);
			grid.setColor(isoColorList[i]);

			MC[1]= ArrayVertices[i+1];
			IsoValue[1] = isoValueList[i+1];
			grid.set(MC[1],IsoValue[1]);
			grid.setColor(isoColorList[i+1]);

			MC[2]= ArrayVertices[i+2];
			IsoValue[2] = isoValueList[i+2];
			grid.set(MC[2],IsoValue[2]);
			grid.setColor(isoColorList[i+2]);

			MC[3]= ArrayVertices[i+3];
			IsoValue[3] = isoValueList[i+3];
			grid.set(MC[3],IsoValue[3]);
			grid.setColor(isoColorList[i+3]);

			MC[4]= ArrayVertices[i+4];
			IsoValue[4] = isoValueList[i+4];
			grid.set(MC[4],IsoValue[4]);
			grid.setColor(isoColorList[i+4]);

			MC[5]= ArrayVertices[i+5];
			IsoValue[5] = isoValueList[i+5];
			grid.set(MC[5],IsoValue[5]);
			grid.setColor(isoColorList[i+5]);

			MC[6]= ArrayVertices[i+6];
			IsoValue[6] = isoValueList[i+6];
			grid.set(MC[6],IsoValue[6]);
			grid.setColor(isoColorList[i+6]);

			MC[7]= ArrayVertices[i+7];
			IsoValue[7] = isoValueList[i+7];
			grid.set(MC[7],IsoValue[7]);
			grid.setColor(isoColorList[i+7]);

			grid.divColor();

			stat = MarchingCube(grid, threshold,triangles,polyCount,polyconnet,polyconnetCount,trianglesColor);
		//}
	}

	MFnMesh meta;
	MDGModifier dgModifier;

	MObject newpoly = meta.create(triangles.length(),polyCount.length(),triangles,polyCount,polyconnet,MObject::kNullObj,&stat);
	//meta.createColorSetWithName("myAwesomeColorSet");
	
	meta.setVertexColors(trianglesColor,polyconnet,&dgModifier,MFnMesh::MColorRepresentation::kRGB );

	if(!stat)
	stat.perror("Unble to create mesh");

	meta.updateSurface();

	MString cmd( "sets -e -fe initialShadingGroup " );
	cmd += meta.name()+ "; ";;

    cmd += ( "polyMergeVertex  -d 0.01 -am 1 -ch 1 " );
	cmd += meta.name() + "; ";

	//cmd += ( "polyQuad  -a 30 -kgb 1 -ktb 1 -khe 1 -ws 1 -ch 1 " );
	//cmd += meta.name() + "; ";

	cmd += ( "select -r " );
	cmd += meta.name();
	MGlobal::executeCommand( cmd );


	// Since this class is derived off of MPxCommand, you can use the 
	// inherited methods to return values and set error messages
	//


	setResult( "MayaMetaball command executed!\n" );

	end = clock();
	double result = (double)(end - begin) / CLOCKS_PER_SEC;
	
	MString txt;
	txt = "\n MC time : ";
	txt += result;
	txt+= "\n";
	MGlobal::displayInfo(txt);


	return redoIt();

}


MStatus MayaMetaPlugin::undoIt()
{
return dgMod.undoIt();
}

MStatus MayaMetaPlugin::redoIt()
{
return dgMod.doIt();
}



MStatus initializePlugin( MObject obj )
{
	MFnPlugin plugin( obj, "S.H Park", "1.0", "2014" );
	

	MStatus stat;
	stat = plugin.registerCommand( "MayaMetaPlugin", MayaMetaPlugin::creator,MayaMetaPlugin::newSyntax );
	if ( !stat )
		stat.perror( "registerCommand failed");

	return stat;
}

MStatus uninitializePlugin( MObject obj )
{
	MFnPlugin plugin( obj );

	MStatus	stat;
	stat = plugin.deregisterCommand( "MayaMetaPlugin" );
	if ( !stat )
		stat.perror( "deregisterCommand failed" );

	return stat;
}
