################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/cluster/cluster_main.c \
../src/cluster/cluster_main2.c 

OBJS += \
./src/cluster/cluster_main.o \
./src/cluster/cluster_main2.o 

C_DEPS += \
./src/cluster/cluster_main.d \
./src/cluster/cluster_main2.d 


# Each subdirectory must supply rules for building sources it contributes
src/cluster/%.o: ../src/cluster/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: K1dp C compiler (Andey)'
	k1-gcc -O0 -g3 -Wall -mcore=k1dp  -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


