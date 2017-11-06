#!/usr/bin/env xmipp_python
# -------------------------------------------------------------------------------
# Protocol header automatically generated for program: xmipp_showj
# {begin_of_header}


#----------------------------------------------------------------------------------------------------
# {section}{has_question}
#----------------------------------------------------------------------------------------------------

# Show
show__K_0002_Show = True

#  --input   arglist
"""
Input files to show
"""
_K_0003_P_input_A_arglist = ""

#  --memory   mem
"""
Memory ammount for JVM
"""
_K_0004_P_memory_A_mem = "2g"

# {list_combo}(,image,gallery,metadata,rotspectra) --mode   mode_value
"""
List of params
"""
_K_0005_P_mode_A_mode_value = "image"

# --poll
"""
Keeps checking for changes on input files  (for image mode only!)
"""
_K_0006_P_poll = False

#  --render   label
"""
Activates images rendering (for metadata mode only)
you can pass which label to render, by default the first one that can be visualized
"""
_K_0007_P_render_A_label = "first"

#  --rows   rows
"""
number of rows in table
"""
_K_0008_P_rows_A_rows = ""

#  --columns   columns
"""
number of columns in table
"""
_K_0009_P_columns_A_columns = ""

#  --zoom   zoom
"""
zoom for images
"""
_K_0010_P_zoom_A_zoom = ""

# {list_combo}(,z,y,x,z_pos,y_pos,x_pos) --view   axis
"""
Viewer position (for volumes only)
"""
_K_0011_P_view_A_axis = "z"

# --dont_apply_geo
"""
Does not read geometrical information(for metadata only)
"""
_K_0012_P_dont_apply_geo = False

# --dont_wrap
"""
Does not wrap (for metadata only)
"""
_K_0013_P_dont_wrap = False

# --debug
"""
debug
"""
_K_0014_P_debug = False

# --mask_toolbar
"""
Open mask toolbar (only valid for images)
"""
_K_0015_P_mask_toolbar = False

# {hidden} Program name
"""This is the name of the program to be executed, dont change this!!!"""
ProgramName = "xmipp_showj"

# {hidden} Show expert options
"""If True, expert options will be displayed """
ShowExpertOptions = False

# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

import sys
from protocol_program import *
from protlib_gui import launchProgramGUI

if __name__ == '__main__':
    launchProgramGUI(sys.argv[0])
