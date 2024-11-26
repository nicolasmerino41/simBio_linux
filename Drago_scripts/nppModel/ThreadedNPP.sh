#!/bin/bash
###################################################################################################################
#############################################  DESCRIPCION ########################################################
###################################################################################################################
#    Este script manda a ejecutar un Script de Julia situada en el directorio indicado                                                                                                                     
#                                                                                                                 
###################################################################################################################
#################################  CONFIGURACION DE LA SOLICITUD AL CLUSTER #######################################
###################################################################################################################
#SBATCH -p generic                	# Particion elegida para el trabajo(cola)
#SBATCH -N 1                        # Numero de nodods solicitados
#SBATCH -n 48     	            # Numero de cores(CPUs)
#SBATCH --mem=188G                 	# memoria total solicitada
#SBATCH -t 0-9:55:00         	    # Duracion solicitada para el trabajo (D-HH:MM:SS)
#SBATCH --job-name=Threaded            # Nombre del trabajo

#SBATCH -o Trabajo-%j.out             # fichero de salida estandart
#SBATCH -e Error-%j.err       # fichero de salida estandart
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
directorio=/lustre/home/mncn/nmerino/nppModel            # Directorio donde se encuentran los Scripts de bash y Julia
cd $directorio                                	      # Nos ubicamos en el directorio de trabajo 
rm -r resultados_threaded
mkdir resultados_threaded                                      # carpeta para almacenar resultados
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

# Fecha/Hora de inicio del trabajo
echo 'Trabajo iniciado en fecha:'
date

# Ejecuta el Script de julia
julia -t 48 threaded_for_npp.jl
#julia -t $ncores run_for_distributed_all_in_one.jl
# Fecha/Hora de fin del trabajo
echo 'Trabajo finalizado en fecha:'
date
