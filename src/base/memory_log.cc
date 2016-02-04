#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "memory_log.h"

//static members
std::auto_ptr<MMU> MMU::_instance;
int MMU::_vmsize =0;
int MMU::_vmpeak =0;
int MMU::_vmrss =0;
int MMU::_vmhwm =0;


MMU * MMU::instance()
{
  if( !_instance.get())
  {
    _instance.reset( new MMU);
  }
  return _instance.get();
}


#ifdef LINUX

#include <unistd.h> //getpid

void MMU::measure()
{
  char buf[2048];
  pid_t mypid = ::getpid();
  snprintf(buf, 2048, "/proc/%d/status", mypid);
  _statistic(buf);
}


int MMU::_statistic(char *pidstatus)
{
	char *line;
	char *vmsize_str;
	char *vmpeak_str;
	char *vmrss_str;
	char *vmhwm_str;
	
	size_t len;
	int    peak;
	
	FILE *f;

	vmsize_str = NULL;
	vmpeak_str = NULL;
	vmrss_str = NULL;
	vmhwm_str = NULL;
	line = (char *)malloc(128);
	len = 128;
	
	f = fopen(pidstatus, "r");
	if (!f) return 1;
	
	/* Read memory size data from /proc/pid/status */
	while (!vmsize_str || !vmpeak_str || !vmrss_str || !vmhwm_str)
	{
		if (getline(&line, &len, f) == -1)
		{
			/* Some of the information isn't there, die */
			return 1;
		}
		
		/* Find VmPeak */
		if (!strncmp(line, "VmPeak:", 7))
		{
			vmpeak_str = strdup(&line[7]);
		}
		
		/* Find VmSize */
		else if (!strncmp(line, "VmSize:", 7))
		{
			vmsize_str = strdup(&line[7]);
		}
		
		/* Find VmRSS */
		else if (!strncmp(line, "VmRSS:", 6))
		{
			vmrss_str = strdup(&line[7]);
		}
		/* Find VmHWM */
		else if (!strncmp(line, "VmHWM:", 6))
		{
			vmhwm_str = strdup(&line[7]);
		}
	}
	free(line);
	
	fclose(f);

	/* Get rid of " kB\n"*/
	len = strlen(vmsize_str);
	vmsize_str[len - 4] = 0;
	len = strlen(vmpeak_str);
	vmpeak_str[len - 4] = 0;
	len = strlen(vmrss_str);
	vmrss_str[len - 4] = 0;
	len = strlen(vmhwm_str);
	vmhwm_str[len - 4] = 0;
	
	/* Output results to stderr */
	//fprintf(stderr, "%s\t%s\t%s\t%s\n", vmsize_str, vmpeak_str, vmrss_str, vmhwm_str);
	
	_vmpeak = atoi(vmpeak_str);
	_vmsize = atoi(vmsize_str);
	_vmrss  = atoi(vmrss_str);
	_vmhwm  = atoi(vmhwm_str);
	
	free(vmpeak_str);
	free(vmsize_str);
	free(vmrss_str);
	free(vmhwm_str);
	
	/* Success */
	return 0;
}
#endif


