#include "flexi.h"

!==================================================================================================================================
!> Unit test 'DeviceMemTools'
!> Test the routines from module: 'DeviceMem'.
!> This file holds the host side code. Device side code in DeviceMemTools.cu
!> @author Spencer Starr - 06.2024
!==================================================================================================================================
MODULE TestDeviceMemTools
USE ISO_C_BINDING
USE MOD_Device
USE MOD_DeviceMem
USE MOD_DeviceManage
IMPLICIT NONE

! Interfaces to C++ methods needed for tests on device side
! C++ implementation of these methods can be found in ./unitTests/DeviceMemTools.cu
INTERFACE
    INTEGER FUNCTION CheckDeviceVarValues(realArrKey, intArrKey) BIND(C, name="CheckDeviceVarValues")
        USE ISO_C_BINDING, ONLY: C_INT
        IMPLICIT NONE
        INTEGER(C_INT),VALUE :: realArrKey
        INTEGER(C_INT),VALUE :: intArrKey
    END FUNCTION CheckDeviceVarValues
END INTERFACE

INTERFACE
    SUBROUTINE ModifyValsOnDevice(realArrKey, intArrKey) BIND(C, name="ModifyValsOnDevice")
        USE ISO_C_BINDING, ONLY: C_INT
        IMPLICIT NONE
        INTEGER(C_INT),VALUE :: realArrKey
        INTEGER(C_INT),VALUE :: intArrKey
    END SUBROUTINE ModifyValsOnDevice
END INTERFACE

INTEGER(C_INT):: d_TestDbleArr, d_TestIntArr

REAL :: TestDbleArr(1000000)
INTEGER :: TestIntArr(1000000)
TYPE(DeviceMemoryState) :: InitialDeviceMemState, DeviceMemStateAfterAlloc, FinalDeviceMemState

CONTAINS

    !==================================================================================================================================
    !> @brief Initialize the test arrays for following test methods
    !==================================================================================================================================
    SUBROUTINE InitTest
    !=========================================================================================

        TestDbleArr = 11.0
        TestIntArr = 51

        ! Use set device call to create the CUDA context now.
        call SetDevice(0)

    END SUBROUTINE InitTest

    !==================================================================================================================================
    !> @brief Test that the device memory allocation methods are working correctly. Implicitly tests method for checking device memory state
    !==================================================================================================================================
    SUBROUTINE CheckDeviceMemAllocation
    ! LOCAL VARIABLES
    INTEGER :: ActualAllocatedDeviceMem
    TYPE(DeviceMemoryState) :: TempMemState
    !=========================================================================================

        ! Check device memory
        CALL CheckDeviceMemoryState(InitialDeviceMemState%FreeMemory, InitialDeviceMemState%TotalMemory)

        ! Allocate double float memory on the device
        CALL AllocateDeviceMemory(d_TestDbleArr, SIZE_C_DOUBLE, 1000000)

        ! Check again to make sure the proper amount of memory was allocated
        CALL CheckDeviceMemoryState(DeviceMemStateAfterAlloc%FreeMemory, DeviceMemStateAfterAlloc%TotalMemory)
        ActualAllocatedDeviceMem = InitialDeviceMemState%FreeMemory - DeviceMemStateAfterAlloc%FreeMemory
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
        IF (ActualAllocatedDeviceMem /= 0) THEN
            PRINT *, "[ERROR - CheckDeviceMemAllocation] Free memory volume on device should be unchanged after hybrid double array allocation."
#else
        IF (ActualAllocatedDeviceMem <= 0) THEN
            PRINT *, "[ERROR - CheckDeviceMemAllocation] Free memory volume on device should be lower after double array allocation."
#endif
            PRINT *, "Device initial memory state     : ", InitialDeviceMemState%FreeMemory
            PRINT *, "Device memory state after 1st alloc : ", DeviceMemStateAfterAlloc%FreeMemory
            ERROR STOP
        END IF

        ! Allocate 4-byte int memory on the device
        CALL AllocateDeviceMemory(d_TestIntArr, SIZE_C_INT, 1000000)
        
        ! Check again to make sure the proper amount of memory was allocated
        TempMemState = DeviceMemStateAfterAlloc
        CALL CheckDeviceMemoryState(DeviceMemStateAfterAlloc%FreeMemory, DeviceMemStateAfterAlloc%TotalMemory)
        ActualAllocatedDeviceMem = TempMemState%FreeMemory - DeviceMemStateAfterAlloc%FreeMemory
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
        IF (ActualAllocatedDeviceMem /= 0) THEN
           PRINT *, "[ERROR - CheckDeviceMemAllocation] Free memory volume on device should be unchanged after hybrid integer array allocation."
#else
        IF (ActualAllocatedDeviceMem <= 0) THEN
           PRINT *, "[ERROR - CheckDeviceMemAllocation] Free memory volume on device should be lower after integer array allocation."
#endif              
           PRINT *, "Device initial memory state     : ", InitialDeviceMemState%FreeMemory
           PRINT *, "Device memory state after 2nd alloc : ", DeviceMemStateAfterAlloc%FreeMemory
           ERROR STOP
        END IF
        
    END SUBROUTINE CheckDeviceMemAllocation


    !==================================================================================================================================
    !> @brief Test methods for copy to device, as well as implicitly test macro for kernel invocation
    !==================================================================================================================================
    SUBROUTINE CheckCopyToDevice
    ! LOCAL VARIABLES
    INTEGER :: deviceErr
    !=========================================================================================

        ! Copy data to the device
        CALL CopyToDevice(d_TestDbleArr, TestDbleArr, 1000000)
        CALL CopyToDevice(d_TestIntArr, TestIntArr, 1000000)

        ! Check to make sure the values on the device are changed
        deviceErr = CheckDeviceVarValues(d_TestDbleArr, d_TestIntArr)

        IF (deviceErr /= 0) THEN
            IF (deviceErr == 1) THEN
                ERROR STOP "[ERROR - CheckCopyToDevice] Real values were not copied to the device correctly"
            ELSE IF (deviceErr == 2) THEN
                ERROR STOP "[ERROR - CheckCopyToDevice] Integer values were not copied to the device correctly"
            ELSE IF (deviceErr == 3) THEN
                ERROR STOP "[ERROR - CheckCopyToDevice] Real and Integer values were not copied to the device correctly"
            ELSE IF (deviceErr == 999) THEN
                ERROR STOP "[ERROR - CheckCopyToDevice] There was an error copying the error code from the device."
            ELSE
                print *, "[ERROR - CheckCopyToDevice] Device error code is: ", deviceErr
                ERROR STOP "[ERROR - CheckCopyToDevice] Kernel returned an unhandled error code"
            END IF
        END IF

    END SUBROUTINE CheckCopyToDevice

    !==================================================================================================================================
    !> @brief Test methods for copy from device
    !==================================================================================================================================
    SUBROUTINE CheckCopyFromDevice
    ! LOCAL VARIABLES
    INTEGER :: deviceErr
    !=========================================================================================

        ! Modify the data on the device
        CALL ModifyValsOnDevice(d_TestDbleArr, d_TestIntArr)

        ! Copy the data back to the host
        CALL CopyFromDevice(TestDbleArr, d_TestDbleArr, 1000000)
        CALL CopyFromDevice(TestIntArr, d_TestIntArr, 1000000)

        ! Check to make sure that the changed values were correctly copied
        IF (TestDbleArr(27) /= 25.0) THEN
           print *, TestDbleArr(27)
            ERROR STOP "Value for double array was not copied back to the host correctly."
        END IF

        IF (TestIntArr(81) /= 45.0) THEN
            ERROR STOP "Value for integer array was not copied back to the host correctly."
        END IF

    END SUBROUTINE CheckCopyFromDevice

    !==================================================================================================================================
    !> @brief Test method for deallocating device memory
    !==================================================================================================================================
    SUBROUTINE CheckDeviceMemFree
    ! LOCAL VARIABLES
    INTEGER :: ActualAllocatedDeviceMem
    !=========================================================================================

        CALL CheckDeviceMemoryState(DeviceMemStateAfterAlloc%FreeMemory, DeviceMemStateAfterAlloc%TotalMemory)

        ! Free the memory on the device
        CALL FreeDeviceMemory(d_TestDbleArr)
        CALL FreeDeviceMemory(d_TestIntArr)

        ! Check device memory again to see that memory was correctly freed
        ! Current memory state on the GPU still stored in DeviceMemStateAfterAlloc
        CALL CheckDeviceMemoryState(FinalDeviceMemState%FreeMemory, FinalDeviceMemState%TotalMemory)
#if defined(USE_HYBRID) && (USE_ACCEL == ACCEL_HIP)
        IF (FinalDeviceMemState%FreeMemory /= DeviceMemStateAfterAlloc%FreeMemory) THEN
           PRINT *, "[ERROR - CheckDeviceMemFree] Final device memory volume after memory free call should be unchanged on hybrid chips."
#else
        IF (FinalDeviceMemState%FreeMemory <= DeviceMemStateAfterAlloc%FreeMemory) THEN
           PRINT *, "[ERROR - CheckDeviceMemFree] Final device memory volume after memory free should be greater than allocated volume."
#endif       
           PRINT *, "Device memory state after alloc : ", DeviceMemStateAfterAlloc%FreeMemory
           PRINT *, "Device final memory state       : ", FinalDeviceMemState%FreeMemory
           ERROR STOP
        END IF
    END SUBROUTINE CheckDeviceMemFree

END MODULE TestDeviceMemTools


!==================================================================================================================================
!> @brief Test driver program
!==================================================================================================================================
PROGRAM DeviceMemTools_TestDriver
USE TestDeviceMemTools
USE MOD_DeviceManage
IMPLICIT NONE

    CALL InitTest
    CALL CheckDeviceMemAllocation
    CALL CheckCopyToDevice
    CALL CheckCopyFromDevice
    CALL CheckDeviceMemFree

END PROGRAM DeviceMemTools_TestDriver
