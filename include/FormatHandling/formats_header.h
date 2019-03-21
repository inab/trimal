// CMake generated code
// Do not manually modify this file

// This file includes all States found, and also defines the ReadWriteMS constructor.
// To be able to be automatically recognized, the new state should:
//        1.- Have the same Class Name as File Name (without the extension)
//        2.- Name must end with '_state'
//        3.- Be placed on ReadWriteMS folder

#include "FormatHandling/clustal_state.h"
#include "FormatHandling/fasta_m10_state.h"
#include "FormatHandling/fasta_state.h"
#include "FormatHandling/htmlreport_state.h"
#include "FormatHandling/mega_interleaved_state.h"
#include "FormatHandling/mega_sequential_state.h"
#include "FormatHandling/nexus_m10_state.h"
#include "FormatHandling/nexus_state.h"
#include "FormatHandling/phylip32_m10_state.h"
#include "FormatHandling/phylip32_state.h"
#include "FormatHandling/phylip40_m10_state.h"
#include "FormatHandling/phylip40_state.h"
#include "FormatHandling/phylip_paml_m10_state.h"
#include "FormatHandling/phylip_paml_state.h"
#include "FormatHandling/pir_state.h"

#include "FormatHandling/FormatManager.h"

namespace FormatHandling {
FormatManager::FormatManager()
{
	addState(new clustal_state(this));
	addState(new fasta_m10_state(this));
	addState(new fasta_state(this));
	addState(new htmlreport_state(this));
	addState(new mega_interleaved_state(this));
	addState(new mega_sequential_state(this));
	addState(new nexus_m10_state(this));
	addState(new nexus_state(this));
	addState(new phylip32_m10_state(this));
	addState(new phylip32_state(this));
	addState(new phylip40_m10_state(this));
	addState(new phylip40_state(this));
	addState(new phylip_paml_m10_state(this));
	addState(new phylip_paml_state(this));
	addState(new pir_state(this));
};
}