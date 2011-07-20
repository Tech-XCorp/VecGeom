#include "UVoxelFinder.hh"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include "UUtils.hh"
#include "UMultiUnion.hh"

using namespace std;

//______________________________________________________________________________   
UVoxelFinder::UVoxelFinder() : fMultiUnion(0),
                               fBoxes(0),
                               fBoundaries(0),
                               fXSortedBoundaries(0),
                               fXNumBound(0),
                               fYSortedBoundaries(0),
                               fYNumBound(0),
                               fZSortedBoundaries(0),
                               fZNumBound(0),
                               fNumNodesSliceX(0),
                               fNumNodesSliceY(0),
                               fNumNodesSliceZ(0),
                               fMemoryX(0),
                               fMemoryY(0),
                               fMemoryZ(0),      
                               fNx(0),
                               fNy(0),
                               fNz(0)
                               {}

//______________________________________________________________________________   
UVoxelFinder::UVoxelFinder(UMultiUnion* multiUnion) : fMultiUnion(multiUnion),
                                                       fBoxes(0),
                                                       fBoundaries(0),
                                                       fXSortedBoundaries(0),
                                                       fXNumBound(0),
                                                       fYSortedBoundaries(0),
                                                       fYNumBound(0),
                                                       fZSortedBoundaries(0),
                                                       fZNumBound(0),
                                                       fNumNodesSliceX(0),
                                                       fNumNodesSliceY(0),
                                                       fNumNodesSliceZ(0),
                                                       fMemoryX(0),
                                                       fMemoryY(0),
                                                       fMemoryZ(0),      
                                                       fNx(0),
                                                       fNy(0),
                                                       fNz(0)
                                                       {}

//______________________________________________________________________________   
UVoxelFinder::~UVoxelFinder()
{
   if(fMultiUnion) delete fMultiUnion;
   if(fBoxes) delete fBoxes;
   if(fBoundaries) delete fBoundaries;
   if(fXSortedBoundaries) delete fXSortedBoundaries;
   if(fYSortedBoundaries) delete fYSortedBoundaries;
   if(fZSortedBoundaries) delete fZSortedBoundaries;
   if(fNumNodesSliceX) delete fNumNodesSliceX;
   if(fNumNodesSliceY) delete fNumNodesSliceY;
   if(fNumNodesSliceZ) delete fNumNodesSliceZ;   
   if(fMemoryX) delete fMemoryX;
   if(fMemoryY) delete fMemoryY;
   if(fMemoryZ) delete fMemoryZ;   
}

//______________________________________________________________________________   
void UVoxelFinder::BuildVoxelLimits()
{
// "BuildVoxelLimits"'s aim is to store the coordinates of the origin as well as
// the half lengths related to the bounding box of each node.
// These quantities are stored in the array "fBoxes" (6 different values per
// node.
   int iIndex = 0;
   int carNodes = fMultiUnion->GetNumNodes(); // Number of nodes in "fMultiUnion"
   UVector3 tempPointConv,tempPoint;                                                  
   if(!carNodes) return; // If no node, no need to carry on
   int fNBoxes = 6*carNodes; // 6 different quantities are stored per node
   if(fBoxes) delete fBoxes;
   fBoxes = new double[fNBoxes]; // Array which will store the half lengths
                                 // related to a particular node, but also
                                 // the coordinates of its origin                                 
      
   /*const*/ VUSolid *tempSolid = NULL;
   /*const*/ UTransform3D *tempTransform = NULL;
   
   double* min = new double[3];
   double* max = new double[3];
   
   for(iIndex = 0 ; iIndex < carNodes ; iIndex++)   
   {
      tempSolid = fMultiUnion->GetSolid(iIndex);
      tempTransform = fMultiUnion->GetTransform(iIndex);

      tempSolid->Extent(min, max);
      /*UUtils::*/TransformLimits(min, max, tempTransform);                        
      
      // Storage of the positions:
      fBoxes[6*iIndex]   = 0.5*(max[0]-min[0]); // dX
      fBoxes[6*iIndex+1] = 0.5*(max[1]-min[1]); // dY
      fBoxes[6*iIndex+2] = 0.5*(max[2]-min[2]); // dZ
      fBoxes[6*iIndex+3] = tempTransform->fTr[0]; // Ox
      fBoxes[6*iIndex+4] = tempTransform->fTr[1]; // Oy
      fBoxes[6*iIndex+5] = tempTransform->fTr[2]; // Oz
   }   
}

//______________________________________________________________________________   
void UVoxelFinder::DisplayVoxelLimits()
{
// "DisplayVoxelLimits" displays the dX, dY, dZ, oX, oY and oZ for each node
   int carNodes = fMultiUnion->GetNumNodes();
   int iIndex = 0;
   
   for(iIndex = 0 ; iIndex < carNodes ; iIndex++)
   {
      cout << "    -> Node " << iIndex+1 << ":\n       * dX = " << fBoxes[6*iIndex] << " ; * dY = " << fBoxes[6*iIndex+1] << " ; * dZ = " << fBoxes[6*iIndex+2] << "\n       * oX = " << fBoxes[6*iIndex+3] << " ; * oY = " << fBoxes[6*iIndex+4] << " ; * oZ = " << fBoxes[6*iIndex+5] << "\n";
   }
}

//______________________________________________________________________________   
void UVoxelFinder::CreateBoundaries()
{
// "CreateBoundaries"'s aim is to determine the slices induced by the bounding boxes,
// along each axis. The created boundaries are stored in the array "fBoundaries"
   int iIndex = 0;
   int carNodes = fMultiUnion->GetNumNodes(); // Number of nodes in structure of
                                            // "UMultiUnion" type
   
   // Determination of the 3D extent of the "UMultIUnion" structure:
   double minima_multi[3], maxima_multi[3];
   fMultiUnion->Extent(minima_multi,maxima_multi);
   
   // Determination of the boundries along x, y and z axis:
   fBoundaries = new double[6*carNodes];
   
   for(iIndex = 0 ; iIndex < carNodes ; iIndex++)   
   {
      // For each node, the boundaries are created by using the array "fBoxes"
      // built in method "BuildVoxelLimits":

      // x boundaries
      fBoundaries[2*iIndex] = fBoxes[6*iIndex+3]-fBoxes[6*iIndex];
      fBoundaries[2*iIndex+1] = fBoxes[6*iIndex+3]+fBoxes[6*iIndex];
      // y boundaries
      fBoundaries[2*iIndex+2*carNodes] = fBoxes[6*iIndex+4]-fBoxes[6*iIndex+1];
      fBoundaries[2*iIndex+2*carNodes+1] = fBoxes[6*iIndex+4]+fBoxes[6*iIndex+1];
      // z boundaries
      fBoundaries[2*iIndex+4*carNodes] = fBoxes[6*iIndex+5]-fBoxes[6*iIndex+2];
      fBoundaries[2*iIndex+4*carNodes+1] = fBoxes[6*iIndex+5]+fBoxes[6*iIndex+2];
   }
}

//______________________________________________________________________________   
void UVoxelFinder::SortBoundaries()
{
// "SortBoundaries" orders the boundaries along each axis (increasing order)
// and also does not take into account redundant boundaries, ie if two boundaries
// are separated by a distance strictly inferior to "UUtils::kTolerance".
// The sorted boundaries are respectively stored in:
//              * fXSortedBoundaries
//              * fYSortedBoundaries
//              * fZSortedBoundaries
// In addition, the number of elements contained in the three latter arrays are
// precised thanks to variables: fXNumBound, fYNumBound and fZNumBound.
   int iIndex = 0;
   int carNodes = fMultiUnion->GetNumNodes();
   int *indexSortedBound = new int[2*carNodes];
   double *tempBoundaries = new double[2*carNodes];
   
   // x axis:
   // -------
   int number_boundaries = 0; // Total number of boundaries after treatment:
   UUtils::Sort(2*carNodes,&fBoundaries[0],&indexSortedBound[0],false);
   
   for(iIndex = 0 ; iIndex < 2*carNodes ; iIndex++)
   {
      if(number_boundaries == 0)
      {
         tempBoundaries[number_boundaries] = fBoundaries[indexSortedBound[iIndex]];
         number_boundaries++;
         continue;
      }

      // If two successive boundaries are too close from each other, only the first one is considered      
      if(UUtils::Abs(tempBoundaries[number_boundaries-1]-fBoundaries[indexSortedBound[iIndex]]) > UUtils::kTolerance)
      {
         tempBoundaries[number_boundaries] = fBoundaries[indexSortedBound[iIndex]];
         number_boundaries++;         
      }
   }
   
   if(fXSortedBoundaries) delete [] fXSortedBoundaries;
   fXSortedBoundaries = new double[number_boundaries];
   memcpy(fXSortedBoundaries,&tempBoundaries[0],number_boundaries*sizeof(double));
   fXNumBound = number_boundaries;
   
   // y axis:
   // -------
   number_boundaries = 0; // Total number of boundaries after treatment:
   UUtils::Sort(2*carNodes,&fBoundaries[2*carNodes],&indexSortedBound[0],false);
   
   for(iIndex = 0 ; iIndex < 2*carNodes ; iIndex++)
   {
      if(number_boundaries == 0)
      {
         tempBoundaries[number_boundaries] = fBoundaries[2*carNodes+indexSortedBound[iIndex]];
         number_boundaries++;
         continue;
      }
      
      if(UUtils::Abs(tempBoundaries[number_boundaries-1]-fBoundaries[2*carNodes+indexSortedBound[iIndex]]) > UUtils::kTolerance)
      {
         tempBoundaries[number_boundaries] = fBoundaries[2*carNodes+indexSortedBound[iIndex]];
         number_boundaries++;         
      }
   }
   
   if(fYSortedBoundaries) delete fYSortedBoundaries;
   fYSortedBoundaries = new double[number_boundaries];
   memcpy(fYSortedBoundaries,&tempBoundaries[0],number_boundaries*sizeof(double));
   fYNumBound = number_boundaries;
   
   // z axis:
   // -------
   number_boundaries = 0; // Total number of boundaries after treatment:
   UUtils::Sort(2*carNodes,&fBoundaries[4*carNodes],&indexSortedBound[0],false);
   
   for(iIndex = 0 ; iIndex < 2*carNodes ; iIndex++)
   {
      if(number_boundaries == 0)
      {
         tempBoundaries[number_boundaries] = fBoundaries[4*carNodes+indexSortedBound[iIndex]];
         number_boundaries++;
         continue;
      }
      
      if(UUtils::Abs(tempBoundaries[number_boundaries-1]-fBoundaries[4*carNodes+indexSortedBound[iIndex]]) > UUtils::kTolerance)
      {
         tempBoundaries[number_boundaries] = fBoundaries[4*carNodes+indexSortedBound[iIndex]];
         number_boundaries++;         
      }
   }
   
   if(fZSortedBoundaries) delete fZSortedBoundaries;
   fZSortedBoundaries = new double[number_boundaries];
   memcpy(fZSortedBoundaries,&tempBoundaries[0],number_boundaries*sizeof(double));
   fZNumBound = number_boundaries;         
}

//______________________________________________________________________________   
void UVoxelFinder::DisplayBoundaries()
{
// Prints the positions of the boundaries of the slices on the three axis:
   int iIndex = 0;
   
   cout << " * X axis:\n    | ";
   for(iIndex = 0 ; iIndex < fXNumBound ; iIndex++)
   {
      cout << fXSortedBoundaries[iIndex] << " ";
      if(iIndex != fXNumBound-1) cout << "-> ";
   }
   cout << "|\n    Number of boundaries: " << fXNumBound << "\n";
   
   cout << " * Y axis:\n    | ";
   for(iIndex = 0 ; iIndex < fYNumBound ; iIndex++)
   {
      cout << fYSortedBoundaries[iIndex] << " ";
      if(iIndex != fYNumBound-1) cout << "-> ";
   }
   cout << "|\n    Number of boundaries: " << fYNumBound << "\n"; 
      
   cout << " * Z axis:\n    | ";
   for(iIndex = 0 ; iIndex < fZNumBound ; iIndex++)
   {
      cout << fZSortedBoundaries[iIndex] << " ";
      if(iIndex != fZNumBound-1) printf("-> ");
   }      
   cout << "|\n    Number of boundaries: " << fZNumBound << "\n"; 
}

//______________________________________________________________________________   
void UVoxelFinder::BuildListNodes()
{
// "BuildListNodes" stores in the arrays "fMemoryX", "fMemoryY" and "fMemoryZ" the
// solids present in each slice aong the three axis.
// In array "fNumNodesSliceX", "fNumNodesSliceY" and "fNumNodesSliceZ" are stored the number of
// solids in the considered slice.
   const double localTolerance = 10E-10;
   int iIndex, jIndex = 0;
   int carNodes = fMultiUnion->GetNumNodes();   
   int nmaxslices = 2*carNodes+1;
   int nperslice = 1+(carNodes-1)/(8*sizeof(char));
   int current = 0;
   
   // Number of slices per axis:
   int xNumslices = fXNumBound-1;
   int yNumslices = fYNumBound-1;
   int zNumslices = fZNumBound-1;

   double xbmin, xbmax, ybmin, ybmax, zbmin, zbmax = 0;
   
   // Memory:
   char *storage = new char[nmaxslices*nperslice];
   memset(storage,0,(nmaxslices*nperslice)*sizeof(char));
   char *bits;
   int number_byte;
   int position;
   
   // Loop on x slices:
   fNumNodesSliceX = new int[xNumslices];
   
   for(iIndex = 0 ; iIndex < xNumslices ; iIndex++)
   {
      bits = &storage[current];
   
      // Loop on the nodes:
      for(jIndex = 0 ; jIndex < carNodes ; jIndex++)
      {
         // Determination of the minimum and maximum position along x
         // of the bounding boxe of each node:
         xbmin = fBoxes[6*jIndex+3]-fBoxes[6*jIndex];
         xbmax = fBoxes[6*jIndex+3]+fBoxes[6*jIndex];
         
         // Is the considered node inside slice?
         if((xbmin - fXSortedBoundaries[iIndex+1]) > -localTolerance) continue;
         if((xbmax - fXSortedBoundaries[iIndex]) < localTolerance) continue;
         
         // The considered node is contained in the current slice:
         fNumNodesSliceX[iIndex]++;
         
         // Storage of the number of the node in the array:
         number_byte = jIndex/8;
         position = jIndex%8;
         bits[number_byte] |= (1 << position);                  
      }
      
//      if(fNumNodesSliceX[iIndex]>0) current += nperslice;
      current += nperslice;
   }
   fNx = current;
   fMemoryX = new char[fNx];
   memcpy(fMemoryX,storage,current*sizeof(char));
   
   // Loop on y slices:
   fNumNodesSliceY = new int[yNumslices];
   current = 0;
   memset(storage,0,(nmaxslices*nperslice)*sizeof(char));
   
   for(iIndex = 0 ; iIndex < yNumslices ; iIndex++)
   {
      bits = &storage[current];
   
      // Loop on the nodes:
      for(jIndex = 0 ; jIndex < carNodes ; jIndex++)
      {
         // Determination of the minimum and maximum position along x
         // of the bounding boxe of each node:
         ybmin = fBoxes[6*jIndex+4]-fBoxes[6*jIndex+1];
         ybmax = fBoxes[6*jIndex+4]+fBoxes[6*jIndex+1];
         
         // Is the considered node inside slice?
         if((ybmin - fYSortedBoundaries[iIndex+1]) > -localTolerance) continue;
         if((ybmax - fYSortedBoundaries[iIndex]) < localTolerance) continue;
         
         // The considered node is contained in the current slice:
         fNumNodesSliceY[iIndex]++;
         
         // Storage of the number of the node in the array:
         number_byte = jIndex/8;
         position = jIndex%8;
         bits[number_byte] |= (1 << position);                  
      }
      
//      if(fNumNodesSliceY[iIndex]>0) current += nperslice;
      current += nperslice;
   }
   fNy = current;
   fMemoryY = new char[fNy];
   memcpy(fMemoryY,storage,current*sizeof(char));    
   
   // Loop on z slices:
   fNumNodesSliceZ = new int[zNumslices];
   current = 0;
   memset(storage,0,(nmaxslices*nperslice)*sizeof(char));
   
   for(iIndex = 0 ; iIndex < zNumslices ; iIndex++)
   {
      bits = &storage[current];
   
      // Loop on the nodes:
      for(jIndex = 0 ; jIndex < carNodes ; jIndex++)
      {
         // Determination of the minimum and maximum position along x
         // of the bounding boxe of each node:
         zbmin = fBoxes[6*jIndex+5]-fBoxes[6*jIndex+2];
         zbmax = fBoxes[6*jIndex+5]+fBoxes[6*jIndex+2];
         
         // Is the considered node inside slice?
         if((zbmin - fZSortedBoundaries[iIndex+1]) > -localTolerance) continue;
         if((zbmax - fZSortedBoundaries[iIndex]) < localTolerance) continue;
         
         // The considered node is contained in the current slice:
         fNumNodesSliceZ[iIndex]++;
         
         // Storage of the number of the node in the array:
         number_byte = jIndex/8;
         position = jIndex%8;
         bits[number_byte] |= (1 << position);                  
      }
      
//      if(fNumNodesSliceZ[iIndex]>0) current += nperslice;
      current += nperslice;
   }
   fNz = current;
   fMemoryZ = new char[fNz];
   memcpy(fMemoryZ,storage,current*sizeof(char));      
}

//______________________________________________________________________________   
void UVoxelFinder::GetCandidatesAsString(const char* mask, std::string &result)
{
// Decode the candidates in mask as string.
   static char buffer[30];
   int carNodes = fMultiUnion->GetNumNodes();   
   char mask1 = 1;
   result = "";
   
   for(int icand=0; icand<carNodes; icand++)
   {
      int byte = icand/8;
      int bit = icand - 8*byte;
   
      if (mask1<<bit & mask[byte])
      {
         sprintf(buffer, "%d ", icand+1); 
         result += buffer;
      }   
   }
}

//______________________________________________________________________________   
void UVoxelFinder::DisplayListNodes()
{
// Prints which solids are present in the slices previously elaborated.
   int iIndex, jIndex = 0;
   int carNodes = fMultiUnion->GetNumNodes();   
   int nperslice = 1+(carNodes-1)/(8*sizeof(char));
   string result = "";
   
   cout << " * X axis:\n";
   for(iIndex = 0 , jIndex = 0 ; iIndex < fXNumBound-1 ; iIndex++ , jIndex += nperslice)
   {
      cout << "    Slice #" << iIndex+1 << ": [" << fXSortedBoundaries[iIndex] << " ; " << fXSortedBoundaries[iIndex+1] << "] -> ";
      
      GetCandidatesAsString(&fMemoryX[jIndex], result);
      cout << "[ " << result.c_str() << "]  \n";
   }
   
   cout << " * Y axis:\n";
   for(iIndex = 0 , jIndex = 0 ; iIndex < fYNumBound-1 ; iIndex++ , jIndex += nperslice)
   {
      cout << "    Slice #" << iIndex+1 << ": [" << fYSortedBoundaries[iIndex] << " ; " << fYSortedBoundaries[iIndex+1] << "] -> ";
      
      GetCandidatesAsString(&fMemoryY[jIndex], result);
      cout << "[ " << result.c_str() << "]  \n";
   }
   
   cout << " * Z axis:\n";
   for(iIndex = 0 , jIndex = 0 ; iIndex < fZNumBound-1 ; iIndex++ , jIndex += nperslice)
   {
      cout << "    Slice #" << iIndex+1 << ": [" << fZSortedBoundaries[iIndex] << " ; " << fZSortedBoundaries[iIndex+1] << "] -> ";
      
      GetCandidatesAsString(&fMemoryZ[jIndex], result);
      cout << "[ " << result.c_str() << "]  \n";
   }    
}

//______________________________________________________________________________   
void UVoxelFinder::Voxelize()
{
   BuildVoxelLimits();
   CreateBoundaries();    
   SortBoundaries(); 
   BuildListNodes();  
}

//______________________________________________________________________________
void UVoxelFinder::TransformLimits(double *min, double *max, UTransform3D *transformation)
{
   int jIndex, kIndex = 0;
   double *vertices = new double[24]; 
   UVector3 tempPointConv,tempPoint;
   double currentX, currentY, currentZ = 0;
   double miniX, miniY, miniZ, maxiX, maxiY, maxiZ;          

   // Detemination of the vertices thanks to the extension of each solid:
      // 1st vertice:
   vertices[ 0] = min[0]; vertices[ 1] = min[1]; vertices[ 2] = min[2];
      // 2nd vertice:
   vertices[ 3] = min[0]; vertices[ 4] = max[1]; vertices[ 5] = min[2];   
      // etc.:
   vertices[ 6] = max[0]; vertices[ 7] = max[1]; vertices[ 8] = min[2];
   vertices[ 9] = max[0]; vertices[10] = min[1]; vertices[11] = min[2];
   vertices[12] = min[0]; vertices[13] = min[1]; vertices[14] = max[2];
   vertices[15] = min[0]; vertices[16] = max[1]; vertices[17] = max[2];
   vertices[18] = max[0]; vertices[19] = max[1]; vertices[20] = max[2];
   vertices[21] = max[0]; vertices[22] = min[1]; vertices[23] = max[2];   
   
   // Loop on th vertices
   for(jIndex = 0 ; jIndex < 8 ; jIndex++)
   {
      kIndex = 3*jIndex;
      tempPoint.Set(vertices[kIndex],vertices[kIndex+1],vertices[kIndex+2]);
      // From local frame to the gobal one:
      tempPointConv = transformation->GlobalPoint(tempPoint);
      
      // Current positions on the three axis:         
      currentX = tempPointConv.x;
      currentY = tempPointConv.y;
      currentZ = tempPointConv.z;
      
      // Initialization of extrema:
      if(jIndex == 0)
      {
         miniX = maxiX = currentX;
         miniY = maxiY = currentY;
         miniZ = maxiZ = currentZ;
         continue;                  
      }
         
      // If need be, replacement of the min & max values:
      if(currentX > maxiX)
         maxiX = currentX;
      if(currentX < miniX)
         miniX = currentX;

      if(currentY > maxiY)
         maxiY = currentY;
      if(currentY < miniY)
         miniY = currentY;  

      if(currentZ > maxiZ)
         maxiZ = currentZ;
      if(currentZ < miniZ)
         miniZ = currentZ;                             
   }
   // Recopy of the extrema in the passed pointers:
   min[0] = miniX;
   min[1] = miniY;
   min[2] = miniZ;
   max[0] = maxiX;
   max[1] = maxiY;
   max[2] = maxiZ; 
} 
