#!/bin/bash
###################################################################################################################
#############################################  DESCRIPCION ########################################################
###################################################################################################################
#    Este script manda a ejecutar un Script de Julia situada en el directorio indicado                                                                                                                     
#                                                                                                                 
###################################################################################################################
#################################  CONFIGURACION DE LA SOLICITUD AL CLUSTER #######################################
###################################################################################################################
#SBATCH -p long         	# Particion elegida para el trabajo(cola)
#SBATCH -N 8                        # Numero de nodods solicitados
#SBATCH --ntasks-per-node=2     	            # Numero de cores(CPUs)
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000Mb                	# memoria total solicitada
#SBATCH -t 3-16:00:00         	    # Duracion solicitada para el trabajo (D-HH:MM:SS)
#SBATCH --job-name=Dist_julia            # Nombre del trabajo

#SBATCH -o TrabajoDistribuido-%j.out             # fichero de salida estandart
#SBATCH -e ErrorDistribuido-%j.err       # fichero de salida estandart

###################################################################################################################
################################# DESPLIEGUE DEL SOFTWARE A UTILIZAR  #############################################
###################################################################################################################

# Limpiamos los posibles modulos que tengamos cargados previamente para evitar problemas
module purge

# Cargamos la rama de software deseada
module load rama0.3

# Cargamos los modulos necesarios para ejecutar el trabajo
module load Julia/1.10.5-linux-x86_64

###################################################################################################################
################################################ FICHEROS #########################################################
###################################################################################################################
directorio=/lustre/home/mncn/nmerino/simBio            # Directorio donde se encuentran los Scripts de bash y Julia
cd $directorio                                	      # Nos ubicamos en el directorio de trabajo 
rm -r resultados_distribuidos
mkdir resultados_distribuidos                                      # carpeta para almacenar resultados
cd /lustre/home/mncn/nmerino/simBio/resultados_distribuidos
mkdir outputs
cd $directorio

###################################################################################################################
########################################## COMANDO A EJECUTAR #####################################################
###################################################################################################################
# Informacion sobre el trabajo en Drago
echo ""
echo "#################################################################"
echo "DATOS DE RECURSOS DEL TRABAJO EN DRAGO"
echo "Nombre del Trabajo: $SLURM_JOB_NAME"
echo "Numero de Trabajo: $SLURM_JOB_ID"
echo "Cola de Trabajo: $SLURM_JOB_PARTITION"
echo "Numero de nodos: $SLURM_NNODES"
echo "Numero de Cores (Tareas): $SLURM_NTASKS "
echo "Memoria por Core: $SLURM_MEM_PER_CPU"
echo "Directorio: $SLURM_SUBMIT_DIR" 
echo "#################################################################"
echo ""

# Variables
ncores=$SLURM_NTASKS

# Fecha/Hora de inicio del trabajo
echo 'Trabajo iniciado en fecha:'
date

# Ejecuta el Script de julia
julia --project -t 1 run_for_distributed.jl
#julia -t $ncores run_for_distributed_all_in_one.jl
# Fecha/Hora de fin del trabajo
echo 'Trabajo finalizado en fecha:'
date
