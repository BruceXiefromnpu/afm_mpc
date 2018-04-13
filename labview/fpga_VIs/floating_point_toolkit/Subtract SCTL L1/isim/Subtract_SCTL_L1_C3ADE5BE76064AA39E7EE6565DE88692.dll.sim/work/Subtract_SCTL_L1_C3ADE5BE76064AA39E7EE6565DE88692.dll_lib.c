/**********************************************************************/
/*   ____  ____                                                       */
/*  /   /\/   /                                                       */
/* /___/  \  /                                                        */
/* \   \   \/                                                       */
/*  \   \        Copyright (c) 2003-2009 Xilinx, Inc.                */
/*  /   /          All Right Reserved.                                 */
/* /---/   /\                                                         */
/* \   \  /  \                                                      */
/*  \___\/\___\                                                    */
/***********************************************************************/

#define XSI_DESIGN_DLL_EXPORT true
#include "xsi.h"

struct XSI_INFO xsi_info;

char *IEEE_P_2592010699;
char *STD_STANDARD;
char *IEEE_P_3499444699;
char *XILINXCORELIB_P_2514345707;
char *IEEE_P_1242562249;
char *STD_TEXTIO;
char *XILINXCORELIB_P_1837083571;
char *XILINXCORELIB_P_0652569459;


XSI_DESIGN_DLLESPEC xsiHandle xsi_open(p_xsi_setup_info setup_info)
{
    int argc = 3;
    char* argv[4] = {"Subtract_SCTL_L1_C3ADE5BE76064AA39E7EE6565DE88692.dll", "-nolog", "-sysgen"};
    if (setup_info->logFileName) {
        argc = 4;
        argv[1] = "-log";
        argv[2] = setup_info->logFileName;
        argv[3] = "-sysgen";
    }
    xsi_init_design(argc, argv);
    xsi_register_info(&xsi_info);

    xsi_register_min_prec_unit(-12);
    ieee_p_2592010699_init();
    ieee_p_3499444699_init();
    xilinxcorelib_p_2514345707_init();
    ieee_p_1242562249_init();
    std_textio_init();
    xilinxcorelib_p_1837083571_init();
    xilinxcorelib_p_0652569459_init();
    xilinxcorelib_a_1139420180_3212880686_init();
    xilinxcorelib_a_1710931151_3212880686_init();
    work_a_1902283046_0862517153_init();


    xsi_register_tops("work_a_1902283046_0862517153");

    IEEE_P_2592010699 = xsi_get_engine_memory("ieee_p_2592010699");
    xsi_register_ieee_std_logic_1164(IEEE_P_2592010699);
    STD_STANDARD = xsi_get_engine_memory("std_standard");
    IEEE_P_3499444699 = xsi_get_engine_memory("ieee_p_3499444699");
    XILINXCORELIB_P_2514345707 = xsi_get_engine_memory("xilinxcorelib_p_2514345707");
    IEEE_P_1242562249 = xsi_get_engine_memory("ieee_p_1242562249");
    STD_TEXTIO = xsi_get_engine_memory("std_textio");
    XILINXCORELIB_P_1837083571 = xsi_get_engine_memory("xilinxcorelib_p_1837083571");
    XILINXCORELIB_P_0652569459 = xsi_get_engine_memory("xilinxcorelib_p_0652569459");

    return xsi_run_sysgen_simulation(argc, argv);

}
