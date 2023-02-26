
#include "./headers/output.hpp"
#include <stdio.h>

void OutputVTK(int &nout, LBM &lb)
{
    int Nx = lb.getNx();
    int Ny = lb.getNy();
    int Nz = lb.getNz();
	int		i,j,k;
	char	filename[128];
	FILE	*fp;
	unsigned int array_size;
	unsigned long int offset=0;
	// short  num16; // Int16 2byte
	float  val32; // Float32 4byte

	sprintf(filename,"./field%06d.vtr",nout);
	fp=fopen(filename,"wb");

	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
	fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",Nx,Ny,Nz);
	fprintf(fp,"  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",Nx,Ny,Nz);
	fprintf(fp,"    <PointData>\n");
	fprintf(fp,"    </PointData>\n");
	fprintf(fp,"    <CellData>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Density\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx)*(Ny)*(Nz);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx)*(Ny)*(Nz);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CellType\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx)*(Ny)*(Nz);
	fprintf(fp,"    </CellData>\n");
	fprintf(fp,"    <Coordinates>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateX\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+1);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateY\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Ny+1);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateZ\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nz+1);
	fprintf(fp,"    </Coordinates>\n");
	fprintf(fp,"  </Piece>\n");
	fprintf(fp,"  </RectilinearGrid>\n");
	fprintf(fp,"  <AppendedData encoding=\"raw\">");
	fprintf(fp,"_");

    // Density (cell)
	array_size=1*4*(Nx)*(Ny)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				val32=(float)lb.fluid1[i][j][k].rho; fwrite(&val32,sizeof(float),1,fp);
			}
		}
	}
    // Velocity (cell)
	array_size=3*4*(Nx)*(Ny)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				val32=(float)lb.fluid1[i][j][k].u; fwrite(&val32,sizeof(float),1,fp);
				val32=(float)lb.fluid1[i][j][k].v; fwrite(&val32,sizeof(float),1,fp);
				val32=0.0; fwrite(&val32,sizeof(float),1,fp);
			}
		}
	}

    // CellType (cell)
	array_size=1*4*(Nx)*(Ny)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				val32=(int)lb.fluid1[i][j][k].type; fwrite(&val32,sizeof(int),1,fp);
			}
		}
	}
	// Coordinates (vertices)
	array_size=1*4*(Nx+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(i=0;i<Nx+1;i++){ val32=(float)(i*dx); fwrite(&val32,sizeof(float),1,fp); }

	array_size=1*4*(Ny+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(j=0;j<Ny+1;j++){ val32=(float)(j*dx); fwrite(&val32,sizeof(float),1,fp); }

	array_size=1*4*(Nz+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz+1;k++){ val32=(float)(k*dx); fwrite(&val32,sizeof(float),1,fp); }

	fprintf(fp,"  </AppendedData>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
}


void printLogo()
{
	std::cout << R"(
             _       ____   __  __ 
            | |     |  _ \ |  \/  |
            | |     | |_) || \  / |
            | |     |  _ < | |\/| |
            | |____ | |_) || |  | |
            |______||____/ |_|  |_|
    )";
    
    std::cout << std :: endl << "      - Flow Diagnostics Laboratory ITB - " << std::endl << std::endl;
}