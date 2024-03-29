add_subdirectory(tests)

set(TESTSOURCES ${TESTSOURCES} PARENT_SCOPE)

set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/additionalOutputInputFile.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/calcWaterTemp.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/calib_basins.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/calib_param.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/calibration.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/clcl.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/climate.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/climateYear.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/configFile.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/daily.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/geo.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/globals.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/glacierYear.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/gw_frac.cpp
  	    ${CMAKE_CURRENT_SOURCE_DIR}/initializeWGHM.cpp
  	    ${CMAKE_CURRENT_SOURCE_DIR}/integrateWGHM.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/json11.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lai.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/land.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/permafrost.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/random.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rout_prepare.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/routing.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/s_max.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/snowInElevationFile.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/timestring.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/upstream_stations.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/watergap.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/wghmStateFile.cpp
        PARENT_SCOPE
        )
set(HEADERS
        ${HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/additionalOutputInputFile.h
      	${CMAKE_CURRENT_SOURCE_DIR}/calcWaterTemp.h
        ${CMAKE_CURRENT_SOURCE_DIR}/calib_basins.h
        ${CMAKE_CURRENT_SOURCE_DIR}/calib_param.h
        ${CMAKE_CURRENT_SOURCE_DIR}/calibration.h
        ${CMAKE_CURRENT_SOURCE_DIR}/clcl.h
        ${CMAKE_CURRENT_SOURCE_DIR}/climate.h
        ${CMAKE_CURRENT_SOURCE_DIR}/climateYear.h
	    ${CMAKE_CURRENT_SOURCE_DIR}/common.h
        ${CMAKE_CURRENT_SOURCE_DIR}/configFile.h
        ${CMAKE_CURRENT_SOURCE_DIR}/conversion.h
        ${CMAKE_CURRENT_SOURCE_DIR}/daily.h
        ${CMAKE_CURRENT_SOURCE_DIR}/def.h
        ${CMAKE_CURRENT_SOURCE_DIR}/exception.h
        ${CMAKE_CURRENT_SOURCE_DIR}/geo.h
	    ${CMAKE_CURRENT_SOURCE_DIR}/globals.h
        ${CMAKE_CURRENT_SOURCE_DIR}/glacierYear.h
        ${CMAKE_CURRENT_SOURCE_DIR}/grid.h
        ${CMAKE_CURRENT_SOURCE_DIR}/gw_frac.h
   	    ${CMAKE_CURRENT_SOURCE_DIR}/initializeWGHM.h
	    ${CMAKE_CURRENT_SOURCE_DIR}/integrateWGHM.h
        ${CMAKE_CURRENT_SOURCE_DIR}/json11.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lai.h
        ${CMAKE_CURRENT_SOURCE_DIR}/land.h
        ${CMAKE_CURRENT_SOURCE_DIR}/option.h
        ${CMAKE_CURRENT_SOURCE_DIR}/permafrost.h
        ${CMAKE_CURRENT_SOURCE_DIR}/random.h
        ${CMAKE_CURRENT_SOURCE_DIR}/rout_prepare.h
        ${CMAKE_CURRENT_SOURCE_DIR}/routing.h
        ${CMAKE_CURRENT_SOURCE_DIR}/routing_checks_global.h
        ${CMAKE_CURRENT_SOURCE_DIR}/s_max.h
        ${CMAKE_CURRENT_SOURCE_DIR}/snowInElevationFile.h
        ${CMAKE_CURRENT_SOURCE_DIR}/stack.h
        ${CMAKE_CURRENT_SOURCE_DIR}/timestring.h
        ${CMAKE_CURRENT_SOURCE_DIR}/upstream_stations.h
        ${CMAKE_CURRENT_SOURCE_DIR}/watergap_checks_global.h
        ${CMAKE_CURRENT_SOURCE_DIR}/wghmStateFile.h
        PARENT_SCOPE
        )
set(SOURCEPDAF
        ${SOURCEPDAF}
        ${CMAKE_CURRENT_SOURCE_DIR}/additionalOutputInputFile.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/calcWaterTemp.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/calib_basins.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/calib_param.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/calibration.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/clcl.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/climate.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/climateYear.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/configFile.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/daily.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/enKF2wghmState.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/extractsub.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/geo.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/globals.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/gw_frac.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/initExtractsub.cpp
  	    ${CMAKE_CURRENT_SOURCE_DIR}/initializeWGHM.cpp
  	    ${CMAKE_CURRENT_SOURCE_DIR}/integrateWGHM.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/json11.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lai.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/land.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/matrix.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/option.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/parameterJsonFile.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/permafrost.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/random.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/rmFiles.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rout_prepare.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/routing.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/s_max.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/snowInElevationFile.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/timestring.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/upstream_stations.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/wghmStateFile.cpp
        PARENT_SCOPE
        )
set(HEADERSPDAF
        ${HEADERSPDAF}
        ${CMAKE_CURRENT_SOURCE_DIR}/additionalOutputInputFile.h
      	${CMAKE_CURRENT_SOURCE_DIR}/calcWaterTemp.h
        ${CMAKE_CURRENT_SOURCE_DIR}/calib_basins.h
        ${CMAKE_CURRENT_SOURCE_DIR}/calib_param.h
        ${CMAKE_CURRENT_SOURCE_DIR}/calibration.h
        ${CMAKE_CURRENT_SOURCE_DIR}/clcl.h
        ${CMAKE_CURRENT_SOURCE_DIR}/climate.h
        ${CMAKE_CURRENT_SOURCE_DIR}/climateYear.h
	    ${CMAKE_CURRENT_SOURCE_DIR}/common.h
        ${CMAKE_CURRENT_SOURCE_DIR}/configFile.h
        ${CMAKE_CURRENT_SOURCE_DIR}/conversion.h
        ${CMAKE_CURRENT_SOURCE_DIR}/daily.h
        ${CMAKE_CURRENT_SOURCE_DIR}/def.h
        ${CMAKE_CURRENT_SOURCE_DIR}/enKF2wghmState.h
        ${CMAKE_CURRENT_SOURCE_DIR}/exception.h
        ${CMAKE_CURRENT_SOURCE_DIR}/extractsub.h
        ${CMAKE_CURRENT_SOURCE_DIR}/geo.h
	    ${CMAKE_CURRENT_SOURCE_DIR}/globals.h
        ${CMAKE_CURRENT_SOURCE_DIR}/grid.h
        ${CMAKE_CURRENT_SOURCE_DIR}/gw_frac.h
        ${CMAKE_CURRENT_SOURCE_DIR}/initExtractsub.h
   	    ${CMAKE_CURRENT_SOURCE_DIR}/initializeWGHM.h
	    ${CMAKE_CURRENT_SOURCE_DIR}/integrateWGHM.h
        ${CMAKE_CURRENT_SOURCE_DIR}/json11.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lai.h
        ${CMAKE_CURRENT_SOURCE_DIR}/land.h
        ${CMAKE_CURRENT_SOURCE_DIR}/matrix.h
        ${CMAKE_CURRENT_SOURCE_DIR}/option.h
	    ${CMAKE_CURRENT_SOURCE_DIR}/parameterJsonFile.h
        ${CMAKE_CURRENT_SOURCE_DIR}/permafrost.h
        ${CMAKE_CURRENT_SOURCE_DIR}/random.h
        ${CMAKE_CURRENT_SOURCE_DIR}/rmFiles.h
        ${CMAKE_CURRENT_SOURCE_DIR}/rout_prepare.h
        ${CMAKE_CURRENT_SOURCE_DIR}/routing_checks_global.h
        ${CMAKE_CURRENT_SOURCE_DIR}/routing.h
        ${CMAKE_CURRENT_SOURCE_DIR}/s_max.h
        ${CMAKE_CURRENT_SOURCE_DIR}/snowInElevationFile.h
        ${CMAKE_CURRENT_SOURCE_DIR}/stack.h
        ${CMAKE_CURRENT_SOURCE_DIR}/timestring.h
        ${CMAKE_CURRENT_SOURCE_DIR}/upstream_stations.h
        ${CMAKE_CURRENT_SOURCE_DIR}/watergap_checks_global.h
        ${CMAKE_CURRENT_SOURCE_DIR}/wghmStateFile.h
        PARENT_SCOPE
        )

