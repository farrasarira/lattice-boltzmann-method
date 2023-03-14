
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

void OutputKeEns(int &step, LBM &lb)
{	
	double e_kinetic = 0;
	double enstro = 0;
	double uu = 0 ;
	double vort = 0;

	for (int i = 0; i < lb.getNx(); i++)
	{
		for (int j = 0; j < lb.getNy(); j++)
		{
			for (int k = 0; k < lb.getNz(); k++)
			{
				if (lb.fluid1[i][j][k].type == TYPE_F)
				{ 
					uu = lb.fluid1[i][j][k].u*lb.fluid1[i][j][k].u + lb.fluid1[i][j][k].v*lb.fluid1[i][j][k].v + lb.fluid1[i][j][k].w*lb.fluid1[i][j][k].w;
					e_kinetic = e_kinetic + lb.fluid1[i][j][k].rho * uu;
					
					int i_a = i + 1;
					int i_b = i - 1;
					int j_a = j + 1;
					int j_b = j - 1;
					int k_a = k + 1;
					int k_b = k - 1;

					if (i_b < 1) i_b = lb.getNx()-2;
					else if(i_a > lb.getNx()-2) i_a = 1;

					if (j_b < 1) j_b = lb.getNy()-2;
					else if(j_a > lb.getNy()-2) j_a = 1;

					if (k_b < 1) k_b = lb.getNz()-2;
					else if(k_a > lb.getNz()-2) k_a = 1;
					

					double uy = (lb.fluid1[i][j_a][k].u - lb.fluid1[i][j_b][k].u) / 2;
					double uz = (lb.fluid1[i][j][k_a].u - lb.fluid1[i][j][k_b].u) / 2;
					
					double vx = (lb.fluid1[i_a][j][k].v - lb.fluid1[i_b][j][k].v) / 2;
					double vz = (lb.fluid1[i][j][k_a].v - lb.fluid1[i][j][k_b].v) / 2;

					double wx = (lb.fluid1[i_a][j][k].w - lb.fluid1[i_b][j][k].w) / 2;
					double wy = (lb.fluid1[i][j_a][k].w - lb.fluid1[i][j_b][k].w) / 2;

					vort = (wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy);
					enstro = enstro + lb.fluid1[i][j][k].rho * vort;
				}
			}
		}
	}
	
	std::cout << step << ",";
	std::cout << e_kinetic << ",";
	std::cout << enstro << std::endl;
	

}

void calcError(int &t,LBM &lb)
{
	double kx = 2.0*M_PI/NX;
	double ky = 2.0*M_PI/NY;

	double ua = 0.0;
	double va = 0.0;

	double sumue2 = 0.0;
	double sumve2 = 0.0;

	double sumua2 = 0.0;
	double sumva2 = 0.0;

	double td = 1.0/(NU*(kx*kx+ky*ky));

	for (int i = 0; i < lb.getNx(); i++)
	{
		for (int j = 0; j < lb.getNy(); j++)
		{
			for (int k = 0; k < lb.getNz(); k++)
			{
				if (lb.fluid1[i][j][k].type == TYPE_F)
				{
					double X = i ;
					double Y = j ;
					ua =-0.04 * cos(kx*X) * sin(ky*Y) *exp(-1.0*t/td);
					va = 0.04 * sin(kx*X) * cos(ky*Y) *exp(-1.0*t/td);

					sumue2 += (lb.fluid1[i][j][k].u - ua) * (lb.fluid1[i][j][k].u - ua);
					sumve2 += (lb.fluid1[i][j][k].v - va) * (lb.fluid1[i][j][k].v - va);
					
					sumua2 += ua * ua;
					sumva2 += va * va;
				}
			}
		}
	}

	std::cout << "u error		= " << sqrt(sumue2/sumua2) << std::endl;
	std::cout << "v error		= " << sqrt(sumve2/sumva2) << std::endl; 

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