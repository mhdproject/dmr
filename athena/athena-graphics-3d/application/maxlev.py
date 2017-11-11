import chombovis
import os
import sys
c=chombovis.latest()

# Use a Statefile
#c.misc.restoreState(file_name = 'gm_7.state')

# Set the screen size
c.misc.setRenderWidgetSize(900,900)

# Colormap
#c.cmap.loadUserColormap( './blue.cmap' )
#c.cmap.loadUserColormap( '/usr/local/ChomboVis-4.0.1/usr/share/ChomboVis/data/reverse.cmap' )
c.cmap.loadUserColormap( '/usr/local/ChomboVis-4.1.3/usr/share/ChomboVis/data/reverse.cmap' )
#c.cmap.loadUserColormap( '/usr/local/ChomboVis-4.1.7/usr/share/ChomboVis/data/reverse.cmap' )
#c.cmap.loadUserColormap( os.environ['CHOMBOVIS_HOME'] + '/share/ChomboVis/data/reverse.cmap' )
c.cmap.setActiveColormap( 'User' )
c.cmap.setCmappedRangeMax(16.0)
c.cmap.setCmappedRangeMin(0.0)
c.cmap.setBackgroundColor(rgb=(1.0,1.0,1.0))
c.cmap.showLegend()
c.grid.setColor(rgb=(0.0,0.0,0.0))



# Levels of Refinement Visible
c.reader.setVisibleLevelMax( c.reader.getNumLevels()-1 )

# Isosurfaces
c.iso.toggleVisibility( 1 )
c.iso.setMin( 0.0 )
c.iso.setMax( 8.0 )
c.iso.setNum( 30 )
c.iso.setUseGhostCells(1)


c.annotation.addNewCaption(  )
c.annotation.setText("Double Mach Reflection Shock Test")
c.annotation.setSize(14)
c.annotation.setColor(0,0,0)
c.annotation.setFont('arial')
c.annotation.setPosition(0.1,0.4)

# Slice Visibility
#for a in 'x','y','z':
#    c.slice.toggleVisibility( on_off=0, axis=a )


c.misc.guiUpdate()

# Print a ppm (portable pixelmap) of the image
#c.misc.hardCopy ( outfile_name = c.cmd_line.getInfile() + '.ppm')

# Exit chombovis
#sys.exit(0)
