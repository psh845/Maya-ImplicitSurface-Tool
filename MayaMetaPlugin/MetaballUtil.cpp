#include "MetaballUtil.h"

MStatus Cube_Array(
	 float &Xarray,
	 float &Yarray,
	 float &Zarray,
	 float &width,
	 MFloatArray &XArrayList,
	 MFloatArray &YArrayList,
	 MFloatArray &ZArrayList
	)
{
	float startX = 0.0;
	float startY = 0.0;
	float startZ = 0.0;

	for(float x=0, startX = -(Xarray)/2.0; startX <Xarray/2; x++,startX += width)
	{
		XArrayList.append(startX);
	}

	for(float y=0, startY = -(Yarray)/2.0; startY<Yarray/2; y++, startY += width)
	{
		YArrayList.append(startY);
	}
	
	for(float z=0, startZ = -(Zarray)/2.0; startZ<Zarray/2; z++, startZ +=width)
	{
				
		ZArrayList.append(startZ);
	}
		
	

	return MS::kSuccess;
}

MStatus Cube_box(
	 MPointArray &ArrayVertices,
	 MFloatArray &XArrayList,
	 MFloatArray &YArrayList,
	 MFloatArray &ZArrayList,
	 MIntArray &isoCheckList
	)
{
	
	int numvertices =0;
	int numfaces = 0; 
	int numfaceConnects = 0;
	
	MPointArray boxVertices;

	//MPointArray boxVertices;
	MIntArray boxPolyCounts;
	MIntArray boxPolyConnects;

	boxVertices.clear();
	boxPolyCounts.clear();
	boxPolyConnects.clear();

	/*for(int x=0; x<XArrayList.length()-1; x++)
	{
		for(int y=0; y<YArrayList.length()-1; y++)
		{
			for(int z=0; z<ZArrayList.length()-1; z++)
			{
				if(gridCheck(isoCheckList,YArrayList.length(),ZArrayList.length(),x,y,z))
				{
					boxVertices.append(MPoint(XArrayList[x],YArrayList[y],ZArrayList[z]));
					boxVertices.append(MPoint(XArrayList[x+1],YArrayList[y],ZArrayList[z]));
					boxVertices.append(MPoint(XArrayList[x+1],YArrayList[y],ZArrayList[z+1]));
					boxVertices.append(MPoint(XArrayList[x],YArrayList[y],ZArrayList[z+1]));

					boxVertices.append(MPoint(XArrayList[x],YArrayList[y+1],ZArrayList[z]));
					boxVertices.append(MPoint(XArrayList[x+1],YArrayList[y+1],ZArrayList[z]));
					boxVertices.append(MPoint(XArrayList[x+1],YArrayList[y+1],ZArrayList[z+1]));
					boxVertices.append(MPoint(XArrayList[x],YArrayList[y+1],ZArrayList[z+1]));

					numvertices += 8;
					numfaces += 6;
					numfaceConnects += 24;
				}

			}
		}
	}*/

	for(int i=0; i <ArrayVertices.length(); i+=8)
	{
		if(gridCheck(isoCheckList,i))
		{
			boxVertices.append(ArrayVertices[i]);
			boxVertices.append(ArrayVertices[i+1]);
			boxVertices.append(ArrayVertices[i+2]);
			boxVertices.append(ArrayVertices[i+3]);
			boxVertices.append(ArrayVertices[i+4]);
			boxVertices.append(ArrayVertices[i+5]);
			boxVertices.append(ArrayVertices[i+6]);
			boxVertices.append(ArrayVertices[i+7]);

			numvertices += 8;
			numfaces += 6;
			numfaceConnects += 24;
		}
	}
	/* append를 쓸때는 setLength 설정을 하지 않아야한다.
	boxPolyCounts.setLength(numfaces);
	for( i=0; i < numfaces; i++)
	{
		boxPolyCounts[i]= 4;
	}

	또는
	 for( i=0; i < numfaces; i++)
	{
		boxPolyCounts.append(4);
	}
	
	둘중 택1
	
	*/

	int i;
	int j;

	/*boxPolyCounts.setLength(numfaces);
	for( i=0; i < numfaces; i++)
	{
		boxPolyCounts[i]= 4;
	}
	
	boxPolyConnects.setLength( numfaceConnects );

	for( i=0, j=0; i < numfaceConnects; i+=24, j+=8)
	{
		boxPolyConnects[i] = j;
		boxPolyConnects[i+1] = (j+3);
		boxPolyConnects[i+2] = (j+2);
		boxPolyConnects[i+3] = (j+1);

		boxPolyConnects[i+4] =(j+4);
		boxPolyConnects[i+5] =(j+5);
		boxPolyConnects[i+6] =(j+6);
		boxPolyConnects[i+7] =(j+7);

		boxPolyConnects[i+8] =j;
		boxPolyConnects[i+9] =(j+1);
		boxPolyConnects[i+10] =(j+5);
		boxPolyConnects[i+11] =(j+4);

		boxPolyConnects[i+12] =j;
		boxPolyConnects[i+13] =(j+4);
		boxPolyConnects[i+14] =(j+7);
		boxPolyConnects[i+15] =(j+3);

		boxPolyConnects[i+16] =(j+2);
		boxPolyConnects[i+17] =(j+3);
		boxPolyConnects[i+18] =(j+7);
		boxPolyConnects[i+19] =(j+6);

		boxPolyConnects[i+20] =(j+1);
		boxPolyConnects[i+21] =(j+2);
		boxPolyConnects[i+22] =(j+6);
		boxPolyConnects[i+23] =(j+5);
	}
	
	*/

	 for( i=0; i < numfaces; i++)
	{
		boxPolyCounts.append(4);
	}

	
	for( i=0, j=0; i < numfaceConnects/24; i++, j+=8)
	{
		boxPolyConnects.append(j);
		boxPolyConnects.append(j+1);
		boxPolyConnects.append(j+2);
		boxPolyConnects.append(j+3);

		boxPolyConnects.append(j+4);
		boxPolyConnects.append(j+7);
		boxPolyConnects.append(j+6);
		boxPolyConnects.append(j+5);

		boxPolyConnects.append(j+1);
		boxPolyConnects.append(j);
		boxPolyConnects.append(j+4);
		boxPolyConnects.append(j+5);

		boxPolyConnects.append(j);
		boxPolyConnects.append(j+3);
		boxPolyConnects.append(j+7);
		boxPolyConnects.append(j+4);

		boxPolyConnects.append(j+3);
		boxPolyConnects.append(j+2);
		boxPolyConnects.append(j+6);
		boxPolyConnects.append(j+7);

		boxPolyConnects.append(j+2);
		boxPolyConnects.append(j+1);
		boxPolyConnects.append(j+5);
		boxPolyConnects.append(j+6);
	}

	

	MStatus stat = MS::kSuccess;
	MFnMesh cubeMeshFn;
	MObject newPoly = cubeMeshFn.create(numvertices,numfaces,boxVertices,boxPolyCounts,boxPolyConnects,MObject::kNullObj, &stat);
	if(!stat)
		stat.perror("Unble to create mesh");



	MString cmd( "polyMergeVertex  -d 0.01 -am 1 -ch 1 " );
	cmd += cubeMeshFn.name() + "; ";
	cmd += ( "select -r " );
	cmd += cubeMeshFn.name();
	MGlobal::executeCommand( cmd );

	
	cubeMeshFn.updateSurface();

	//MString txt;

	//for(i=0; i<numvertices; i++)
	//{

	//	txt +=  MString(" ")+ boxVertices[i].x + ", " + boxVertices[i].y+ ", " + boxVertices[i].z + "\n";

	//}
	//for(i=0; i<numfaces; i++)
	//{

	//	txt += MString(" ")+ boxPolyCounts[i] + "\n";
	//}

	//	for(i=0; i<numfaceConnects; i++)
	//{

	//	txt +=  MString(" ")+ boxPolyConnects[i] + "\n";
	//}

	
	/*MGlobal::displayInfo( txt );*/

	/*MString txt;
	txt += MString(" ") + boxVertices.length() + "\n";
	for(int i=0; i<boxVertices.length(); i++)
		txt +=  MString(" ")+ boxVertices[i].x + ", " + boxVertices[i].y+ ", " + boxVertices[i].z + "\n";

	MGlobal::displayInfo( txt );*/

	return MS::kSuccess;
}

MStatus sBox() //box sample
{
	const int numfaces = 6; 
	const int numvertices =8;
	const int numfaceConnects = 24;

	MPointArray boxVertices;
	boxVertices.append(MPoint(-1.0,-1.0,-1.0));
	boxVertices.append(MPoint(-1.0,-1.0,1.0));
	boxVertices.append(MPoint(1.0,-1.0,1.0));
	boxVertices.append(MPoint(1.0,-1.0,-1.0));

	boxVertices.append(MPoint(-1.0,1.0,-1.0));
	boxVertices.append(MPoint(-1.0,1.0,1.0));
	boxVertices.append(MPoint(1.0,1.0,1.0));
	boxVertices.append(MPoint(1.0,1.0,-1.0));

	MIntArray boxPolyCounts;
	MIntArray boxPolyConnects;

	int i;
	/* append를 쓸때는 setLength 설정을 하지 않아야한다.*/
	boxPolyCounts.setLength(numfaces);
	for( i=0; i < boxPolyCounts.length(); i++)
	{
		boxPolyCounts[i] = 4;
	}

	int face_connects[numfaceConnects] = 
		{
			0,3,2,1, 
			4,5,6,7,
			0,1,5,4,
			0,4,7,3,
			2,3,7,6,
			1,2,6,5
		};
	boxPolyConnects.setLength( numfaceConnects );
	for( i=0; i < numfaceConnects; i++)
	{
		
		boxPolyConnects[i] = face_connects[i];
	}

	MStatus stat;
	MFnMesh cubeMeshFn;
	MObject newPoly = cubeMeshFn.create(numvertices,numfaces,boxVertices,boxPolyCounts,boxPolyConnects,MObject::kNullObj, &stat);
	if(!stat)
		stat.perror("Unble to create mesh");

	cubeMeshFn.updateSurface();

	MString cmd( "sets -e -fe initialShadingGroup " );
	cmd += cubeMeshFn.name();
	MGlobal::executeCommand( cmd );

	return MS::kSuccess;
}


float Metaball(
	MPoint &boxVertice, 
	float &width,
	MPointArray &centerSphere,
	MFloatArray &radius,
	MFloatArray &weight,
	MColorArray &color,
	float &isolevel,
	MColor &c
	)
{
	int button = 1;
	float x = boxVertice.x;
	float y = boxVertice.y;
	float z = boxVertice.z;
	float sum= 0.0;
	float w=0.0;
	float d = 0.0;
	float r =0.0;

	

	for (int n=0; n<radius.length(); n++)
	{
		MVector p = boxVertice;

		d = sqrt((centerSphere[n].x - p.x)*(centerSphere[n].x - p.x) + (centerSphere[n].y - p.y)*(centerSphere[n].y - p.y) + (centerSphere[n].z - p.z)* (centerSphere[n].z - p.z));
		//print(weight.get(n)+","+d+":"+ cx.get(n)+","+cy.get(n)+","+cz.get(n) +"\t");
		r = radius[n];

		switch(button)
		{
		case 0:
			/*meta ball*/
			if (d>=0 && d<=r/3)
			{
				w= weight[n]*(1-((3*(d*d))/(r*r)));
				sum+=w;
				ColorMeta(c,w,color[n],isolevel);
			} else if (d>r/3 && d<=r)
			{
				w= (3*weight[n]/2)*((1-d/r)*(1-d/r));
				sum+=w;
				ColorMeta(c,w,color[n],isolevel);
			} else
			{
				w= 0;
				sum+=w;
				//ColorMeta(c,w,color[n],isolevel);
			}
			//print("metaball"+"\t");

			break;

		case 1:
			/*Soft Objects*/
			if (d>=0 && d<=r)
			{
				w= weight[n]*(1-(4*pow(d, 6)/(9*pow(r, 6)))+(17*pow(d, 4)/(9*pow(r, 4)))-(22*pow(d, 2)/(9*pow(r, 2))));
				sum+=w;
				ColorMeta(c,w,color[n],isolevel);
			} else
			{
				w=0;
				sum+=w;
				//ColorMeta(c,w,color[n],isolevel);
			}
			//print("Soft Objects"+"\t");

			break;

		default:
			;
			//print("no metaball type: 0.metaball    1.Soft Objects");
		}
		//print(x+","+y+":"+int(sum)+"\t");
	}
  //print(x+","+y+","+z+":"+sum+"\n");
  //return d;
  return sum;

}


float Metaball_line(
	MPoint &boxVertice, 
	float &width,
	MPointArray &line_point,
	MFloatArray &line_radius,
	MFloatArray &line_radius2,
	MFloatArray &line_weight
	)
{
	int button = 1;
	float x = boxVertice.x;
	float y = boxVertice.y;
	float z = boxVertice.z;
	float sum= 0.0;
	float d = 0.0;
	float r =0.0;
	float r1 =0.0;
	float r2 =0.0;
	


	for (int n=0; n<line_radius.length(); n++)
	{

		MVector p = boxVertice;
		MVector p0 = line_point[n*2];
		MVector p1 = line_point[n*2+1];

		MVector v = p1;
		v.operator-=(p0);
		MVector t = p;
		t.operator-=(p0);
		float a = t.operator*(v)/v.operator*(v);
		v.operator*=(a);
		v.operator+=(p0);

		
		r1 = line_radius[n];
		r2 = line_radius2[n];

		if (a<=0)
		{
		  d = sqrt((p.x - p0.x)*(p.x - p0.x) + (p.y - p0.y)*(p.y - p0.y) + (p.z - p0.z)*(p.z - p0.z));
		  r=r1;
		} else if (a>0 && a<1)
		{
		  d = sqrt((p.x - v.x)*(p.x - v.x) + (p.y - v.y)*(p.y - v.y) + (p.z - v.z)*(p.z - v.z));

		  r= a*r2 + (1-a)*r1;
		} else
		{
		  d = sqrt((p.x - p1.x)*(p.x - p1.x) + (p.y - p1.y)*(p.y - p1.y) + (p.z - p1.z)*(p.z - p1.z));
		  r =r2;
		}

		switch(button)
		{
		case 0:
			/*meta ball*/
			if (d>=0 && d<=r/3)
			{
				sum+= line_weight[n]*(1-((3*(d*d))/(r*r)));
			} else if (d>r/3 && d<=r)
			{
				sum+= (3*line_weight[n]/2)*((1-d/r)*(1-d/r));
			} else
			{
				sum+= 0;
			}
			//print("metaball"+"\t");

			break;

		case 1:
			/*Soft Objects*/
			if (d>=0 && d<=r)
			{
				sum+= line_weight[n]*(1-(4*pow(d, 6)/(9*pow(r, 6)))+(17*pow(d, 4)/(9*pow(r, 4)))-(22*pow(d, 2)/(9*pow(r, 2))));
			} else
			{
				sum+=0;
			}
			//print("Soft Objects"+"\t");

			break;

		default:
			;
			//print("no metaball type: 0.metaball    1.Soft Objects");
		}
		//print(x+","+y+":"+int(sum)+"\t");
	}
  //print(x+","+y+","+z+":"+sum+"\n");
  //return d;
  return sum;

}


float Metaball_curve4p(MPoint &boxVertice,
	float &width,
	MPointArray &curve_point,
	MFloatArray &curve_radius,
	MFloatArray &curve_radius2,
	MFloatArray &curve_weight
	)
{
	int button = 1;
	float x = boxVertice.x;
	float y = boxVertice.y;
	float z = boxVertice.z;
	float sum= 0.0;
	float d = 0.0;
	float dd = 0.0;
	float r =0.0;
	float r1 =0.0;
	float r2 =0.0;
	float rt =0.0;
	float rt2 = 0.0;


	for (int n=0; n<curve_radius.length(); n++)
	{

		MVector p = boxVertice;
		MVector p0 = curve_point[n*4];
		MVector p1 = curve_point[n*4+1];
		MVector p2 = curve_point[n*4+2];
		MVector p3 = curve_point[n*4+3];
		MVector pf = p0;
		MVector pt = p0;

		r1 = curve_radius[n];
		r2 = curve_radius2[n];

		for(float t =0; t<=1.0; t+=0.05)
		{

			t=floor(t*100)*0.01;
			float b0 = (1 - t) * (1 - t) * (1 - t); 
			float b1 = 3 * t * (1 - t) * (1 - t); 
			float b2 = 3 * t * t * (1 - t); 
			float b3 = t * t * t; 

			pt.x = p0.x*b0 + p1.x*b1 + p2.x*b2 + p3.x*b3;
			pt.y = p0.y*b0 + p1.y*b1 + p2.y*b2 + p3.y*b3;
			pt.z = p0.z*b0 + p1.z*b1 + p2.z*b2 + p3.z*b3;

			MVector v = pt;
			v.operator-=(pf);
			MVector tt = p;
			tt.operator-=(pf);

			float a = tt.operator*(v)/v.operator*(v);
			v.operator*=(a);
			v.operator+=(pf);


			if(a<=0)
			{
				dd = sqrt((p.x - pf.x)*(p.x - pf.x) + (p.y - pf.y)*(p.y - pf.y) + (p.z - pf.z)*(p.z - pf.z));
				rt = r1;
			}
			else if(a>0 && a<1)
			{
				dd = sqrt((p.x - v.x)*(p.x - v.x) + (p.y - v.y)*(p.y - v.y) + (p.z - v.z)*(p.z - v.z));
				rt2 = r1+ (r2-r1)*t;
				rt = a*rt2 + (1-a)*r1;
			}
			else
			{
				dd = sqrt((p.x - pt.x)*(p.x - pt.x) + (p.y - pt.y)*(p.y - pt.y) + (p.z - pt.z)*(p.z - pt.z));
				rt = r1+ (r2-r1)*t;
			}


			if(t==0)
			{
				d=dd;
				r=rt;
			}
			else
			{
				if(dd<d){
					d=dd;
					r=rt;
				}
			}

			if(t<=1.0)
			{
				pf.x=pt.x;
				pf.y=pt.y;
				pf.z=pt.z;
				r1 =rt;
			}

		}

		//r = curve_radius[n];
		switch(button)
		{
		case 0:
			/*meta ball*/
			if (d>=0 && d<=r/3)
			{
				sum+= curve_weight[n]*(1-((3*(d*d))/(r*r)));
			} else if (d>r/3 && d<=r)
			{
				sum+= (3*curve_weight[n]/2)*((1-d/r)*(1-d/r));
			} else
			{
				sum+= 0;
			}
			//print("metaball"+"\t");

			break;

		case 1:
			/*Soft Objects*/
			if (d>=0 && d<=r)
			{
				sum+= curve_weight[n]*(1-(4*pow(d, 6)/(9*pow(r, 6)))+(17*pow(d, 4)/(9*pow(r, 4)))-(22*pow(d, 2)/(9*pow(r, 2))));
			} else
			{
				sum+=0;
			}
			//print("Soft Objects"+"\t");

			break;

		default:
			;
			//print("no metaball type: 0.metaball    1.Soft Objects");
		}
		//print(x+","+y+":"+int(sum)+"\t");
	}
  //print(x+","+y+","+z+":"+sum+"\n");
  //return d;
  return sum;

	
}



float Metaball_curve(MPoint &boxVertice,
	float &width,
	MPointArray &curve_point,
	MFloatArray &curve_radius,
	MFloatArray &curve_radius2,
	MFloatArray &curve_weight,
	MIntArray &curve_npoint
	)
{
	int button = 1;
	float x = boxVertice.x;
	float y = boxVertice.y;
	float z = boxVertice.z;
	float sum= 0.0;
	float d = 0.0;
	float dd = 0.0;
	float r =0.0;
	int tt=0;
	float r1 =0.0;
	float r2 =0.0;
	float rt =0.0;
	float rt2 = 0.0;


	for (int n=0; n<curve_radius.length(); n++)
	{

		MVector p = boxVertice;
		MPointArray pp;
		pp.clear();
		
		for(int i=0; i<curve_npoint[n]; i++)
		{
			pp.append(curve_point[tt+i]);
		}

		int ntring = curve_npoint[n]-1;

		MVector pf=pp[0];
		MVector pt;

		r1 = curve_radius[n];
		r2 = curve_radius2[n];

		for(float t =0; t<=1.0; t+=0.05)
		{
			
			int k,kn,nn,nkn;
			float blend,muk,munk;

			t=floor(t*100)*0.01;

			
			pt.x=0.0;
			pt.y=0.0;
			pt.z=0.0;

			muk = 1;
			munk = pow(1-t,(float)ntring);

			 for ( k=0;k<=ntring;k++)
			{
				nn = ntring;
				kn = k;
				nkn = ntring - k;
				blend = muk * munk;
				muk *= t;
				munk /= (1-t);
				while (nn >= 1) 
				{
					blend *= nn;
					nn--;
					if (kn > 1) {
						blend /= (float)kn;
					kn--;
					}
					if (nkn > 1) {
						blend /= (float)nkn;
						nkn--;
					}
				}
				pt.x += pp[k].x * blend;
				pt.y += pp[k].y * blend;
				pt.z += pp[k].z * blend;
			}
			
			

			MVector v = pt;
			v.operator-=(pf);
			MVector tt = p;
			tt.operator-=(pf);

			float a = tt.operator*(v)/v.operator*(v);
			v.operator*=(a);
			v.operator+=(pf);


			if(a<=0)
			{
				dd = sqrt((p.x - pf.x)*(p.x - pf.x) + (p.y - pf.y)*(p.y - pf.y) + (p.z - pf.z)*(p.z - pf.z));
				rt = r1;
			}
			else if(a>0 && a<1)
			{
				dd = sqrt((p.x - v.x)*(p.x - v.x) + (p.y - v.y)*(p.y - v.y) + (p.z - v.z)*(p.z - v.z));
				rt2 = r1+ (r2-r1)*t;
				rt = a*rt2 + (1-a)*r1;
			}
			else
			{
				dd = sqrt((p.x - pt.x)*(p.x - pt.x) + (p.y - pt.y)*(p.y - pt.y) + (p.z - pt.z)*(p.z - pt.z));
				rt = r1+ (r2-r1)*t;
			}


			if(t==0)
			{
				d=dd;
				r=rt;
			}
			else
			{
				if(dd<d){
					d=dd;
					r=rt;
				}
			}

			if(t<=1.0)
			{
				pf.x=pt.x;
				pf.y=pt.y;
				pf.z=pt.z;
				r1 =rt;
			}

		}

		//r = curve_radius[n];
		switch(button)
		{
		case 0:
			/*meta ball*/
			if (d>=0 && d<=r/3)
			{
				sum+= curve_weight[n]*(1-((3*(d*d))/(r*r)));
			} else if (d>r/3 && d<=r)
			{
				sum+= (3*curve_weight[n]/2)*((1-d/r)*(1-d/r));
			} else
			{
				sum+= 0;
			}
			//print("metaball"+"\t");

			break;

		case 1:
			/*Soft Objects*/
			if (d>=0 && d<=r)
			{
				sum+= curve_weight[n]*(1-(4*pow(d, 6)/(9*pow(r, 6)))+(17*pow(d, 4)/(9*pow(r, 4)))-(22*pow(d, 2)/(9*pow(r, 2))));
			} else
			{
				sum+=0;
			}
			//print("Soft Objects"+"\t");

			break;

		default:
			;
			//print("no metaball type: 0.metaball    1.Soft Objects");
		}
		//print(x+","+y+":"+int(sum)+"\t");
		tt +=ntring+1;
	}
  //print(x+","+y+","+z+":"+sum+"\n");
  //return d;
  return sum;

	
}



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
	)
{
	int button = 1;
	float x = boxVertice.x;
	float y = boxVertice.y;
	float z = boxVertice.z;
	float w=0.0;
	float sum= 0.0;
	float d = 0.0;
	float dd = 0.0;
	float r =0.0;
	int tt=0;
	float r1 =0.0;
	float r2 =0.0;
	float rt =0.0;
	float rt2 = 0.0;


	MDagPath dagPath;
	MObject component;
	MFnNurbsCurve curvFn;
	MItSelectionList iter(selection,MFn::kNurbsCurve);
	//iter.setFilter(MFn::kNurbsCurve);

	MVector p = boxVertice;
	int n=0;
	for(;!iter.isDone();iter.next())
	{

		iter.getDagPath(dagPath);
		curvFn.setObject(dagPath);

		double tStart, tEnd;
		curvFn.getKnotDomain(tStart, tEnd);
		int aa =curvFn.numCVs();

		MPoint pp;

		float t;
		float ts;
		//float tincr =(float)(tEnd -tStart)*0.05;
		//float tincr = 0.05*(tEnd -tStart);
		float tincr = 0.05;
		
		curvFn.getPointAtParam(tStart,pp,MSpace::kWorld);
		MVector pf = pp;
		MVector pt = pp;


		/*curve_radius.append(2);
		curve_weight.append(2);

		curve_radius2.append(3);*/

		r1 = curve_radius[n];
		r2 = curve_radius2[n];

		srand((unsigned int)time(NULL));
		float ra =rand()%11*0.1;
		srand(time(NULL));
		float ra2 =rand()%11*0.1;
		srand(time(NULL));
		float ra3 =rand()%11*0.1;
		
		/*if(n==0)
			curve_color.append(1.0,0.0,0.0);
		else if(n==1)
			curve_color.append(0.0,0.0,1.0);
		else
			curve_color.append(ra,ra2,ra3);*/
		//curve_color.append(1.0,0.0,0.0);

		//for(t =tStart; t<=tEnd; t+=tincr)
		for(t =0.0; t<=1.0; t+=tincr)
		{
			t=floor(t*100.)/100;
			ts=t*(tEnd -tStart);
			curvFn.getPointAtParam(ts,pp,MSpace::kWorld);
			pt = pp;

			MVector v = pt;
			v.operator-=(pf);
			MVector tt = p;
			tt.operator-=(pf);

			float a = (float)((tt.operator*(v))/(v.operator*(v)));
			v.operator*=(a);
			v.operator+=(pf);


			if(a<=0.0)
			{
				dd = sqrt((p.x - pf.x)*(p.x - pf.x) + (p.y - pf.y)*(p.y - pf.y) + (p.z - pf.z)*(p.z - pf.z));
				rt = r1;
			}
			else if(a>0.0 && a<1.0)
			{
				dd = sqrt((p.x - v.x)*(p.x - v.x) + (p.y - v.y)*(p.y - v.y) + (p.z - v.z)*(p.z - v.z));
				rt2 = r1+ (r2-r1)*t;
				rt = a*rt2 + (1-a)*r1;
			}
			else
			{
				dd = sqrt((p.x - pt.x)*(p.x - pt.x) + (p.y - pt.y)*(p.y - pt.y) + (p.z - pt.z)*(p.z - pt.z));
				rt = r1+ (r2-r1)*t;
			}


			if(t==0.0){
				d=dd;
				r=rt;
			}
			else
			{
				if(dd<d){
					d=dd;
					r=rt;
				}
			}

			if(t<=1.0)
			{
				pf.x=pt.x;
				pf.y=pt.y;
				pf.z=pt.z;
				r1 =rt;
			}

		}
		
		//r = curve_radius[n];
		switch(button)
		{
		case 0:
			/*meta ball*/
			if (d>=0 && d<=r/3)
			{
				w= curve_weight[n]*(1-((3*(d*d))/(r*r)));
				sum +=w;
				ColorMeta(c,w,curve_color[n],isolevel);

			} else if (d>r/3 && d<=r)
			{
				 w=  (3*curve_weight[n]/2)*((1-d/r)*(1-d/r));
				 sum+=w;
				 ColorMeta(c,w,curve_color[n],isolevel);
			} else
			{
				w= 0;
				sum+=w;
				//ColorMeta(c,w,curve_color[n],isolevel);
			}
			//print("metaball"+"\t");

			break;

		case 1:
			/*Soft Objects*/
			if (d>=0 && d<=r)
			{
				w= curve_weight[n]*(1-(4*pow(d, 6)/(9*pow(r, 6)))+(17*pow(d, 4)/(9*pow(r, 4)))-(22*pow(d, 2)/(9*pow(r, 2))));
				sum+=w;
				ColorMeta(c,w,curve_color[n],isolevel);
			} else
			{
				w=0;
				sum+= w;
				//ColorMeta(c,w,curve_color[n],isolevel);
			}
			//print("Soft Objects"+"\t");

			break;

		default:
			;
			//print("no metaball type: 0.metaball    1.Soft Objects");
			//w=curve_weight[n]*exp(-r*(d*d));//blooby
		}
		//print(x+","+y+":"+int(sum)+"\t");

		n++;
	}


  //print(x+","+y+","+z+":"+sum+"\n");
  //return d;
  return sum;

}


MStatus MarchingCube(
	GridCell &grid,
	float &isolevel,
	MPointArray &triangles,
	MIntArray &polyCount,
	MIntArray &polyconnet,
	int &polyconnetCount,
	MColorArray &trianglesColor
	)
{
	
	int i,ntriang;
	int cubeindex;
	MPoint vertlist[12];
	MColor Colorlist[12];

	int edgeTable[256]={
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

	int triTable[256][16] =
	{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
	{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
	{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
	{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
	{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
	{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
	{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
	{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
	{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
	{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
	{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
	{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
	{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
	{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
	{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
	{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
	{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
	{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
	{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
	{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
	{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
	{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
	{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
	{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
	{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
	{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
	{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
	{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
	{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
	{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
	{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
	{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
	{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
	{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
	{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
	{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
	{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
	{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
	{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
	{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
	{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
	{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
	{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
	{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
	{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
	{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
	{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
	{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
	{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
	{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
	{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
	{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
	{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
	{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
	{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
	{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
	{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
	{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
	{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
	{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
	{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
	{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
	{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
	{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
	{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
	{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
	{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
	{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
	{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
	{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
	{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
	{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
	{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
	{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
	{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
	{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
	{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
	{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
	{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
	{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
	{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
	{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
	{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
	{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
	{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
	{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
	{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
	{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
	{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
	{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
	{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
	{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
	{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
	{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
	{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
	{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
	{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

   /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
   */
   cubeindex = 0;
  if (grid.val[0] > isolevel) cubeindex |= 1;
  if (grid.val[1] > isolevel) cubeindex |= 2;
  if (grid.val[2] > isolevel) cubeindex |= 4;
  if (grid.val[3] > isolevel) cubeindex |= 8;
  if (grid.val[4] > isolevel) cubeindex |= 16;
  if (grid.val[5] > isolevel) cubeindex |= 32;
  if (grid.val[6] > isolevel) cubeindex |= 64;
  if (grid.val[7] > isolevel) cubeindex |= 128;

   /* Cube is entirely in/out of the surface */
   if (edgeTable[cubeindex] == 0)
      return MS::kSuccess;

   /* Find the vertices where the surface intersects the cube */
   if (edgeTable[cubeindex] & 1){
      vertlist[0] =
	  VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);

	  Colorlist[0] =
		ColorInterp(isolevel, grid.c[0], grid.c[1], grid.val[0], grid.val[1]);
   }
   if (edgeTable[cubeindex] & 2){
      vertlist[1] =
         VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);

	  Colorlist[1] =
		ColorInterp(isolevel, grid.c[1], grid.c[2], grid.val[1], grid.val[2]);
   }
   if (edgeTable[cubeindex] & 4){
      vertlist[2] =
         VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);

	  Colorlist[2] =
		ColorInterp(isolevel, grid.c[2], grid.c[3], grid.val[2], grid.val[3]);
   }
   if (edgeTable[cubeindex] & 8){
      vertlist[3] =
         VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);

	  Colorlist[3] =
		ColorInterp(isolevel, grid.c[3], grid.c[0], grid.val[3], grid.val[0]);
   }
   if (edgeTable[cubeindex] & 16){
      vertlist[4] =
         VertexInterp(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);

	  Colorlist[4] =
		ColorInterp(isolevel, grid.c[4], grid.c[5], grid.val[4], grid.val[5]);
   }
   if (edgeTable[cubeindex] & 32){
      vertlist[5] =
         VertexInterp(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);

	  Colorlist[5] =
		ColorInterp(isolevel, grid.c[5], grid.c[6], grid.val[5], grid.val[6]);
   }
   if (edgeTable[cubeindex] & 64){
      vertlist[6] =
         VertexInterp(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);

	  Colorlist[6] =
		ColorInterp(isolevel, grid.c[6], grid.c[7], grid.val[6], grid.val[7]);
   }
   if (edgeTable[cubeindex] & 128){
      vertlist[7] =
         VertexInterp(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);

	  Colorlist[7] =
		ColorInterp(isolevel, grid.c[7], grid.c[4], grid.val[7], grid.val[4]);
   }
   if (edgeTable[cubeindex] & 256){
      vertlist[8] =
         VertexInterp(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);

	  Colorlist[8] =
		ColorInterp(isolevel, grid.c[0], grid.c[4], grid.val[0], grid.val[4]);
   }
   if (edgeTable[cubeindex] & 512){
      vertlist[9] =
         VertexInterp(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);

	  Colorlist[9] =
		ColorInterp(isolevel, grid.c[1], grid.c[5], grid.val[1], grid.val[5]);
   }
   if (edgeTable[cubeindex] & 1024){
      vertlist[10] =
         VertexInterp(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);

	  Colorlist[10] =
		ColorInterp(isolevel, grid.c[2], grid.c[6], grid.val[2], grid.val[6]);
   }
   if (edgeTable[cubeindex] & 2048){
      vertlist[11] =
         VertexInterp(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);

	  Colorlist[11] =
		ColorInterp(isolevel, grid.c[3], grid.c[7], grid.val[3], grid.val[7]);
   }
	

   /* Create the triangle */
   ntriang = 0;
   MPoint p;
   MColor c;
   MColor cc;
  

   for (i=0;triTable[cubeindex][i]!=-1;i+=3) 
   {
	  
	  p = vertlist[triTable[cubeindex][i]];
	  triangles.append(p);
	  polyconnet.append(polyconnetCount++);
		

	  p =vertlist[triTable[cubeindex][i+1]];
	  triangles.append(p);
	  polyconnet.append(polyconnetCount++);

	  p=vertlist[triTable[cubeindex][i+2]];
      triangles.append(p);
	  polyconnet.append(polyconnetCount++);

	  polyCount.append(3);

	  c= (0.0,0.0,0.0);
	  cc= (0.0,0.0,0.0);
	  c=(Colorlist[triTable[cubeindex][i]]);
	  trianglesColor.append(c);
	  cc.operator+=(c);
	  c=(Colorlist[triTable[cubeindex][i+1]]);
	  trianglesColor.append(c);
	  cc.operator+=(c);
	  c=(Colorlist[triTable[cubeindex][i+2]]);
	  trianglesColor.append(c);
	  cc.operator+=(c);
	  cc.operator/=(3);


	  ntriang++;
   }

 /*  MString txt;
	txt += MString(" ") + triangles.length() + "\n";
	for(int i=0; i<triangles.length(); i++)
		txt +=  MString(" ")+ triangles[i].x + ", " + triangles[i].y+ ", " + triangles[i].z + "\n";

	for(int i=0; i<polyCount.length(); i++)
		txt +=  MString(" ")+ polyCount[i] + "\n";


	for(int i=0; i<polyconnet.length(); i++)
		txt +=  MString(" ")+ polyconnet[i] + "\n";

	MGlobal::displayInfo( txt );*/
   

	return MS::kSuccess;
}


MPoint VertexInterp(float &isolevel,MPoint &p1,MPoint &p2, float &valp1, float &valp2)
{
   float mu;
   MPoint p;

   if (abs(isolevel-valp1) < 0.00001)
      return(p1);
   if (abs(isolevel-valp2) < 0.00001)
      return(p2);
   if (abs(valp1-valp2) < 0.00001)
      return(p1);
   mu = (isolevel - valp1) / (valp2 - valp1);
   p.x = p1.x + mu * (p2.x - p1.x);
   p.y = p1.y + mu * (p2.y - p1.y);
   p.z = p1.z + mu * (p2.z - p1.z);

   return(p);
}

MColor ColorInterp(float isolevel, MColor c1, MColor c2, float valp1, float valp2)
{

  float mu;
  MColor c;

  if (abs(isolevel-valp1) < 0.00001)
    return(c1);
  if (abs(isolevel-valp2) < 0.00001)
    return(c2);
  if (abs(valp1-valp2) < 0.00001)
    return(c1);
  mu = (isolevel - valp1) / (valp2 - valp1);
  c.r = c1.r + mu * (c2.r - c1.r);
  c.g = c1.g + mu * (c2.g - c1.g);
  c.b = c1.b + mu * (c2.b - c1.b);

  return(c);
}


void ColorMeta(MColor &c,float w, MColor color,float isolevel)
{

  float mu;
  if (w<=0)
  {
	c.r += 0;
    c.g += 0;
    c.b+= 0;
  } else if (w<=isolevel)
  {
	mu = (w - 0) / (isolevel - 0);
    c.r+= 0 + mu *(color.r-0);
    c.g+= 0 + mu *(color.g-0);
    c.b+= 0 + mu *(color.b-0);
  } else
  {
    c.r+= color.r;
    c.g+= color.g;
    c.b+= color.b;
  }
  //return c;
}


bool gridCheck(
	MIntArray &isoCheckList,
	int val)
{
  bool  boolean=false;



  if(isoCheckList[val]==1)
    boolean= true;

  if (isoCheckList[val+1]==1)
    boolean= true;
  if (isoCheckList[val+2]==1)
    boolean= true;
  if (isoCheckList[val+3]==1)
    boolean= true;

  if (isoCheckList[val+4]==1)
    boolean= true;
  if (isoCheckList[val+5]==1)
    boolean= true;
  if (isoCheckList[val+6]==1)
    boolean= true;
  if (isoCheckList[val+7]==1)
    boolean= true;

  return boolean;
}

MStatus getData_curve(
	MSelectionList &selection,
	MFloatArray &Curve_radius,
	MFloatArray &Curve_radius2,
	MFloatArray &Curve_weight,
	MColorArray &Curve_color
	)
{
	MDagPath dagPath;
	
	MFnNurbsCurve curvFn;
	MFnNurbsSurface surfaceFn;
	MFnTransform transformFn;
	MString name;
	MItSelectionList iter(selection,MFn::kNurbsCurve);

    MFnDependencyNode nodeFn; 
	MStatus stat = MS::kSuccess;

	//iter.setFilter(MFn::kNurbsCurve);
	
	for (;!iter.isDone();iter.next())
    {
		MObject component;
		iter.getDependNode(component);
		nodeFn.setObject(component);
		
		 
		//MGlobal::displayInfo(MString("11")+ nodeFn.name().asChar() );

		MPlug curve_rad1 =nodeFn.findPlug("curve_rad1",&stat);
		if(stat == MStatus::kSuccess)
		{
		float cr1;
		curve_rad1.getValue(cr1);
		//cr1=curve_rad1.asFloat()
		Curve_radius.append(cr1);
		}

		MPlug curve_rad2 =nodeFn.findPlug("curve_rad2",&stat);
		if(stat == MStatus::kSuccess)
		{
		float cr2;
		curve_rad2.getValue(cr2);
		Curve_radius2.append(cr2);
		}

		MPlug curve_w =nodeFn.findPlug("curve_weight",&stat);
		if(stat == MStatus::kSuccess)
		{
		float cw;
		curve_w.getValue(cw);
		Curve_weight.append(cw);
		}

		MColor cc;
		float cr;
		float cg;
		float cb;
		MPlug curve_cX =nodeFn.findPlug("curve_colorX",&stat);
		if(stat == MStatus::kSuccess)
		{
		
			curve_cX.getValue(cr);
	
			MPlug curve_cY =nodeFn.findPlug("curve_colorY",&stat);
		
		
			curve_cY.getValue(cg);
	
			MPlug curve_cZ =nodeFn.findPlug("curve_colorZ",&stat);
		
		
			curve_cZ.getValue(cb);

			cc.r = cr;
			cc.g = cg;
			cc.b = cb;
			Curve_color.append(cc);
		}
		
		

	}


	return stat;
}

MStatus getData(
	MSelectionList &selection,
	MPointArray &centerSphere,
	MFloatArray &radius,
	MFloatArray &weight,
	MColorArray &color
	
	)
{
	MDagPath dagPath;
	
	
	MFnTransform transformFn;
	MString name;
	MItSelectionList iter(selection,MFn::kNurbsSurface);

    MFnDependencyNode nodeFn; 
	MFnMesh nl;
	MStatus stat = MS::kSuccess;

	//iter.setFilter(MFn::kMesh);
	
	for (;!iter.isDone();iter.next())
    {
		MObject component;
		iter.getDependNode(component);
		nodeFn.setObject(component);
		//nl.setObject(component);
		//nl.f
		
		//MGlobal::displayInfo( MString("22")+nodeFn.name().asChar() );
		MPoint p;
		double tx;
		double ty;
		double tz;

		MPlug sph_transX =nodeFn.findPlug("sph_transX",stat);
		if(stat == MStatus::kSuccess)
		{
				sph_transX.getValue(tx);
		
				MPlug sph_transY =nodeFn.findPlug("sph_transY",stat);
			
				sph_transY.getValue(ty);
	
				MPlug sph_transZ =nodeFn.findPlug("sph_transZ",stat);
			
				sph_transZ.getValue(tz);

				p.x = tx;
				p.y = ty;
				p.z = tz;
				centerSphere.append(p);
		}

		MPlug sph_rad =nodeFn.findPlug("sph_rad",stat);
		if(stat == MStatus::kSuccess)
		{
		float rad;
		sph_rad.getValue(rad);
		radius.append(rad);
		}


		MPlug sph_weight =nodeFn.findPlug("sph_weight",stat);
		{
		float w;
		sph_weight.getValue(w);
		weight.append(w);
		}

		MColor c;
		float cr;
		float cg;
		float cb;
		MPlug sph_colorX =nodeFn.findPlug("sph_colorX",stat);
		if(stat == MStatus::kSuccess)
		{
		
			sph_colorX.getValue(cr);
			MPlug sph_colorY =nodeFn.findPlug("sph_colorY",stat);
		
			sph_colorY.getValue(cg);
			MPlug sph_colorZ =nodeFn.findPlug("sph_colorZ",stat);
		
			sph_colorZ.getValue(cb);
			c.r = cr;
			c.g = cg;
			c.b = cb;
			color.append(c);
		
		}
		

	}

	return stat;
}