# The Channel To Wire mapping needs these libraries
# All of this assumes you have set up the icarusalg package so you can interface to LArSoft services
import ROOT
##import ICARUSservices as services
import galleryUtils

# And here recover the necessary services
def getServices(services):
    detClocks      = services.ServiceManager('DetectorClocks')
    detProperties  = services.ServiceManager('DetectorProperties')
    geometryCore   = services.ServiceManager('Geometry')
    lar_properties = services.ServiceManager('LArProperties')

    return detClocks,detProperties,geometryCore,lar_properties
