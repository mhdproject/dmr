#############################################
# MAKE A SINGLE MOVIE FRAME USING CHOMBOVIS #
#############################################

from vtkpython import *
import chombovis
import anag_utils
import sys
import os

def frame_it(dir, state_file):
#c = chombovis.ChomboVis()
    c = chombovis.latest()

    c.misc.restoreState(file_name = dir + '/' + state_file)\

    c.misc.setRenderWidgetSize(800,800)

#    c.misc.closeAllDialogs()

#    print ' no of levels = ',c.reader.getNumberOfLevels()
#    c.reader.setVisibleLevelMax(c.reader.getNumberOfLevels()-1)
    
#    c.misc.vtkUpdate()
    
    c.misc.hardCopy( outfile_name = c.cmd_line.getInfile() + '.ppm', \
                     mag_factor=2, \
                     format_extension='ppm' )
    
    print ' Wrote File ' + c.cmd_line.getInfile() + '.ppm'

    return

###    c.misc.close()

#############################################################################


state_file = 'gm_6.state'

dir = './'

frame_it(dir, state_file)

