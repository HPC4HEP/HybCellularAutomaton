################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/host/host_main.c 

OBJS += \
./src/host/host_main.o 

C_DEPS += \
./src/host/host_main.d 


# Each subdirectory must supply rules for building sources it contributes
src/host/%.o: ../src/host/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: K1dp C compiler (Andey)'
	k1-gcc -O0 -g3 -Wall -mcore=k1dp  -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


