from opentrons import protocol_api

# metadata
metadata = {
    'protocolName': '3-input AND gates, all states and induction matrix',
    'author': 'Rizki Mardian',
    'description': 'Protocol for preparing inducers stock for 3-input AND gates',
    'apiLevel': '2.0'
}

# protocol run function. the part after the colon lets your editor know
# where to look for autocomplete suggestions
def run(protocol: protocol_api.ProtocolContext):

    r_pipette_name = 'p300_single'
    r_tiprack_slots = ['1', '2', '4', '5']
    r_tiprack_name = 'opentrons_96_tiprack_300ul'
    r_tip_racks = [protocol.load_labware(r_tiprack_name, slot) for slot in r_tiprack_slots]
    r_pipette = protocol.load_instrument(instrument_name = r_pipette_name, mount = 'right', tip_racks = r_tip_racks)

    eppendorf = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', '3')
    plate = protocol.load_labware('thermofisher_96_wellplate_450ul', '6')

    counter = 0
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    cols = (4, 11)
    output_row_idx = 0 #start from the first row
    output_col_idx = cols[0] #start from the third col

    volume = 20
    
    ### column 1 - all induction states ###
    for i in range(8):
        
        induction = format(i, '03b')
                
        r_pipette.pick_up_tip()
        r_pipette.aspirate(volume, eppendorf['A{}'.format(int(induction[0])+5)])
        r_pipette.dispense(volume, plate[rows[i]+'1'])
        r_pipette.drop_tip()

        r_pipette.pick_up_tip()
        r_pipette.aspirate(volume, eppendorf['B{}'.format(int(induction[1])+5)])
        r_pipette.dispense(volume, plate[rows[i]+'1'])
        r_pipette.drop_tip()
        
        r_pipette.pick_up_tip()
        r_pipette.aspirate(volume, eppendorf['C{}'.format(int(induction[2])+5)])
        r_pipette.dispense(volume, plate[rows[i]+'1'])
        r_pipette.drop_tip()

    ### column 4:11 - induction matrix ###
    for a in range(4):
        for b in range(4):
            for c in range(4):
                
                if output_col_idx > cols[1]:
                    output_row_idx += 1
                    output_col_idx = cols[0]
                
                r_pipette.pick_up_tip()
                r_pipette.aspirate(volume, eppendorf['A{}'.format(a+1)])
                r_pipette.dispense(volume, plate[rows[output_row_idx]+str(output_col_idx)])
                r_pipette.drop_tip()

                r_pipette.pick_up_tip()
                r_pipette.aspirate(volume, eppendorf['B{}'.format(b+1)])
                r_pipette.dispense(volume, plate[rows[output_row_idx]+str(output_col_idx)])
                r_pipette.drop_tip()

                r_pipette.pick_up_tip()
                r_pipette.aspirate(volume, eppendorf['C{}'.format(c+1)])
                r_pipette.dispense(volume, plate[rows[output_row_idx]+str(output_col_idx)])
                r_pipette.drop_tip()
                
                output_col_idx += 1
                counter += 1