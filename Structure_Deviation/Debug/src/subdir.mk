################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Main.cpp \
../src/biomolecule.cpp \
../src/coordinationsite.cpp \
../src/geometry.cpp \
../src/ligand.cpp \
../src/report.cpp \
../src/shape.cpp \
../src/standard.cpp \
../src/submolecule.cpp \
../src/tetragonalpl_structure_cd.cpp \
../src/tetragonalpl_structure_fe.cpp \
../src/tetragonalpl_structure_mg.cpp \
../src/tetragonalpl_structure_ni.cpp \
../src/tetragonalplanar.cpp \
../src/trigonalpl_structure_as.cpp \
../src/trigonalpl_structure_ca.cpp \
../src/trigonalpl_structure_cd.cpp \
../src/trigonalpl_structure_co.cpp \
../src/trigonalpl_structure_cu.cpp \
../src/trigonalpl_structure_fe.cpp \
../src/trigonalpl_structure_hg.cpp \
../src/trigonalpl_structure_k.cpp \
../src/trigonalpl_structure_mg.cpp \
../src/trigonalpl_structure_mn.cpp \
../src/trigonalpl_structure_na.cpp \
../src/trigonalpl_structure_ni.cpp \
../src/trigonalpl_structure_pb.cpp \
../src/trigonalpl_structure_zn.cpp \
../src/trigonalplanar.cpp 

OBJS += \
./src/Main.o \
./src/biomolecule.o \
./src/coordinationsite.o \
./src/geometry.o \
./src/ligand.o \
./src/report.o \
./src/shape.o \
./src/standard.o \
./src/submolecule.o \
./src/tetragonalpl_structure_cd.o \
./src/tetragonalpl_structure_fe.o \
./src/tetragonalpl_structure_mg.o \
./src/tetragonalpl_structure_ni.o \
./src/tetragonalplanar.o \
./src/trigonalpl_structure_as.o \
./src/trigonalpl_structure_ca.o \
./src/trigonalpl_structure_cd.o \
./src/trigonalpl_structure_co.o \
./src/trigonalpl_structure_cu.o \
./src/trigonalpl_structure_fe.o \
./src/trigonalpl_structure_hg.o \
./src/trigonalpl_structure_k.o \
./src/trigonalpl_structure_mg.o \
./src/trigonalpl_structure_mn.o \
./src/trigonalpl_structure_na.o \
./src/trigonalpl_structure_ni.o \
./src/trigonalpl_structure_pb.o \
./src/trigonalpl_structure_zn.o \
./src/trigonalplanar.o 

CPP_DEPS += \
./src/Main.d \
./src/biomolecule.d \
./src/coordinationsite.d \
./src/geometry.d \
./src/ligand.d \
./src/report.d \
./src/shape.d \
./src/standard.d \
./src/submolecule.d \
./src/tetragonalpl_structure_cd.d \
./src/tetragonalpl_structure_fe.d \
./src/tetragonalpl_structure_mg.d \
./src/tetragonalpl_structure_ni.d \
./src/tetragonalplanar.d \
./src/trigonalpl_structure_as.d \
./src/trigonalpl_structure_ca.d \
./src/trigonalpl_structure_cd.d \
./src/trigonalpl_structure_co.d \
./src/trigonalpl_structure_cu.d \
./src/trigonalpl_structure_fe.d \
./src/trigonalpl_structure_hg.d \
./src/trigonalpl_structure_k.d \
./src/trigonalpl_structure_mg.d \
./src/trigonalpl_structure_mn.d \
./src/trigonalpl_structure_na.d \
./src/trigonalpl_structure_ni.d \
./src/trigonalpl_structure_pb.d \
./src/trigonalpl_structure_zn.d \
./src/trigonalplanar.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


