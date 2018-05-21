--------------------------------------------------------------------------------
--    This file is owned and controlled by Xilinx and must be used solely     --
--    for design, simulation, implementation and creation of design files     --
--    limited to Xilinx devices or technologies. Use with non-Xilinx          --
--    devices or technologies is expressly prohibited and immediately         --
--    terminates your license.                                                --
--                                                                            --
--    XILINX IS PROVIDING THIS DESIGN, CODE, OR INFORMATION "AS IS" SOLELY    --
--    FOR USE IN DEVELOPING PROGRAMS AND SOLUTIONS FOR XILINX DEVICES.  BY    --
--    PROVIDING THIS DESIGN, CODE, OR INFORMATION AS ONE POSSIBLE             --
--    IMPLEMENTATION OF THIS FEATURE, APPLICATION OR STANDARD, XILINX IS      --
--    MAKING NO REPRESENTATION THAT THIS IMPLEMENTATION IS FREE FROM ANY      --
--    CLAIMS OF INFRINGEMENT, AND YOU ARE RESPONSIBLE FOR OBTAINING ANY       --
--    RIGHTS YOU MAY REQUIRE FOR YOUR IMPLEMENTATION.  XILINX EXPRESSLY       --
--    DISCLAIMS ANY WARRANTY WHATSOEVER WITH RESPECT TO THE ADEQUACY OF THE   --
--    IMPLEMENTATION, INCLUDING BUT NOT LIMITED TO ANY WARRANTIES OR          --
--    REPRESENTATIONS THAT THIS IMPLEMENTATION IS FREE FROM CLAIMS OF         --
--    INFRINGEMENT, IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A   --
--    PARTICULAR PURPOSE.                                                     --
--                                                                            --
--    Xilinx products are not intended for use in life support appliances,    --
--    devices, or systems.  Use in such applications are expressly            --
--    prohibited.                                                             --
--                                                                            --
--    (c) Copyright 1995-2018 Xilinx, Inc.                                    --
--    All rights reserved.                                                    --
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- You must compile the wrapper file Multiply_Accumulator_635C28D3DC43420AADC518895231A915.vhd when simulating
-- the core, Multiply_Accumulator_635C28D3DC43420AADC518895231A915. When compiling the wrapper file, be sure to
-- reference the XilinxCoreLib VHDL simulation library. For detailed
-- instructions, please refer to the "CORE Generator Help".

-- The synthesis directives "translate_off/translate_on" specified
-- below are supported by Xilinx, Mentor Graphics and Synplicity
-- synthesis tools. Ensure they are correct for your synthesis tool(s).

LIBRARY ieee;
USE ieee.std_logic_1164.ALL;
-- synthesis translate_off
LIBRARY XilinxCoreLib;
-- synthesis translate_on
ENTITY Multiply_Accumulator_635C28D3DC43420AADC518895231A915 IS
  PORT (
    clk : IN STD_LOGIC;
    ce : IN STD_LOGIC;
    sclr : IN STD_LOGIC;
    bypass : IN STD_LOGIC;
    a : IN STD_LOGIC_VECTOR(15 DOWNTO 0);
    b : IN STD_LOGIC_VECTOR(15 DOWNTO 0);
    s : OUT STD_LOGIC_VECTOR(31 DOWNTO 0)
  );
END Multiply_Accumulator_635C28D3DC43420AADC518895231A915;

ARCHITECTURE Multiply_Accumulator_635C28D3DC43420AADC518895231A915_a OF Multiply_Accumulator_635C28D3DC43420AADC518895231A915 IS
-- synthesis translate_off
COMPONENT wrapped_Multiply_Accumulator_635C28D3DC43420AADC518895231A915
  PORT (
    clk : IN STD_LOGIC;
    ce : IN STD_LOGIC;
    sclr : IN STD_LOGIC;
    bypass : IN STD_LOGIC;
    a : IN STD_LOGIC_VECTOR(15 DOWNTO 0);
    b : IN STD_LOGIC_VECTOR(15 DOWNTO 0);
    s : OUT STD_LOGIC_VECTOR(31 DOWNTO 0)
  );
END COMPONENT;

-- Configuration specification
  FOR ALL : wrapped_Multiply_Accumulator_635C28D3DC43420AADC518895231A915 USE ENTITY XilinxCoreLib.xbip_multaccum_v2_0(behavioral)
    GENERIC MAP (
      c_a_type => 0,
      c_a_width => 16,
      c_accum_mode => 0,
      c_accum_width => 32,
      c_b_type => 0,
      c_b_width => 16,
      c_bypass_low => 0,
      c_ce_overrides_sclr => 0,
      c_has_bypass => 1,
      c_latency => 1,
      c_out_width => 32,
      c_round_type => 0,
      c_use_dsp48 => 1,
      c_verbosity => 0,
      c_xdevicefamily => "spartan6"
    );
-- synthesis translate_on
BEGIN
-- synthesis translate_off
U0 : wrapped_Multiply_Accumulator_635C28D3DC43420AADC518895231A915
  PORT MAP (
    clk => clk,
    ce => ce,
    sclr => sclr,
    bypass => bypass,
    a => a,
    b => b,
    s => s
  );
-- synthesis translate_on

END Multiply_Accumulator_635C28D3DC43420AADC518895231A915_a;
