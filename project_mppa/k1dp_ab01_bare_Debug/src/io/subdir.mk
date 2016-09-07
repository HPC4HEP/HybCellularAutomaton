################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/io/io_main.c \
../src/io/io_main2.c 

OBJS += \
./src/io/io_main.o \
./src/io/io_main2.o 

C_DEPS += \
./src/io/io_main.d \
./src/io/io_main2.d 


# Each subdirectory must supply rules for building sources it contributes
src/io/%.o: ../src/io/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: K1dp C compiler (Andey)'
	k1-gcc -O0 -g3 -Wall -mcore=k1dp  -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


