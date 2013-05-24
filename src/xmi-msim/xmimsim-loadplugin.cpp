
/*
Copyright (C) 2013 Bruno Golosio and Tom Schoonjans

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "xrmc_loadxmimsim.h"
#include <string>
#include <glib/gprintf.h> 
#include <iostream>
#include "xrmc_exception.h"

using namespace std;

GModule *xrmc_xmimsim = NULL;

int LoadXrmcXmimsimPlugin() {

  if (xrmc_xmimsim != NULL) {
  	cout << "XRMC XMI-MSIM plugin alread loaded." << endl;
  }

  //Initialize xrmc_xmimsim
  XmiCheckXrmcXmimsimPlugin xmi_check_xrmc_xmimsim_plugin;
  gchar *module_path;
  gchar *plugin_dir;
  
  if (!g_module_supported()) {
    throw xrmc_exception("GModule not supported\n");
  }
  //environment variable could also be useful here
  if ((plugin_dir = (gchar *) g_getenv("XRMC_XMIMSIM_MODULE")) == NULL)
    plugin_dir = g_strdup(XRMC_XMIMSIM_LIB);
  
  module_path = g_strdup_printf("%s" G_DIR_SEPARATOR_S "%s.%s", plugin_dir, "xrmc-xmimsim", G_MODULE_SUFFIX);
  xrmc_xmimsim = g_module_open(module_path, (GModuleFlags) 0);
  g_free(module_path);
  
  if (!xrmc_xmimsim) {
    g_fprintf(stderr,"GModule error message: %s\n", g_module_error());
    throw xrmc_exception("Could not load module xrmc-xmimsim\n");
  }
  
  if (!g_module_symbol(xrmc_xmimsim, "xmi_check_xrmc_xmimsim_plugin", (gpointer *) &xmi_check_xrmc_xmimsim_plugin)) {
    g_module_close(xrmc_xmimsim);
    g_fprintf(stderr, "GModule error message: %s\n", g_module_error());
    throw xrmc_exception("Could not get symbol xmi_check_xrmc_xmimsim_plugin from module xrmc-xmimsim\n");
  }
  
  if (xmi_check_xrmc_xmimsim_plugin == NULL) {
    g_module_close(xrmc_xmimsim);
    throw xrmc_exception("Symbol xmi_check_xrmc_xmimsim_plugin from module xrmc-xmimsim is NULL\n");
  }
  
  if (xmi_check_xrmc_xmimsim_plugin() == 1) {
    g_fprintf(stdout,"xrmc_xmimsim plugin is functional\n");
  }
  return 0;
}
