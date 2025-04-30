# Overview of code:
'''
1. define setup
2. take parameter for temperature and RPM
3. pause to add CHX, then transfer to cold plate and mix
To simulate:
pixi run opentrons_simulate -L ./labware protocols/Yeast_OPP_v1.py
'''

from opentrons import protocol_api
from opentrons import types
import math, re

metadata = {'protocolName': 'Yeast_OPP_v1','author': 'Jared Bard','source': 'Protocol Library',}
requirements = {"robotType": "Flex","apiLevel": "2.21",}

def add_parameters(parameters):
    # ======================== RUNTIME PARAMETERS ========================
    parameters.add_bool(
        display_name="Dry Run",
        variable_name="DRYRUN",
        default=False,
        description="Whether to perform a dry run or not.")
    parameters.add_bool(
        display_name="Use Temp",
        variable_name="USETEMP",
        default=True,
        description="Use temperature for heater/shaker")
    parameters.add_int(
        display_name="Temp",
        variable_name="TEMP",
        default=37, minimum = 37, maximum = 95,
        description="Temperature for heater/shaker")
    parameters.add_int(
        display_name="RPM",
        variable_name="RPM",
        default=1300, minimum = 100, maximum = 2500,
        description="Speed of heater/shaker")
    
def run(protocol: protocol_api.ProtocolContext):
    # region ======================== DOWNLOADED PARAMETERS ========================
    DRYRUN              = protocol.params.DRYRUN
    USETEMP         = protocol.params.USETEMP
    TEMP             = protocol.params.TEMP
    RPM         = protocol.params.RPM

    TIP50_APINAME = 'opentrons_flex_96_filtertiprack_50ul'
    TIP1000_APINAME = 'opentrons_flex_96_filtertiprack_1000ul'
    # endregion
    # region =============================== PIPETTE ===============================
    p50 = protocol.load_instrument("flex_8channel_50", "left")
    p1000 = protocol.load_instrument('flex_8channel_1000', 'right')
    p1000_flow_rate_aspirate_default = 600
    p1000_flow_rate_dispense_default = 600
    p1000_flow_rate_blow_out_default = 400
    p50_flow_rate_aspirate_default = 50
    p50_flow_rate_dispense_default = 50
    p50_flow_rate_blow_out_default = 100
    # endregion
    # region =============================== Deck Setup ===============================

    temp_block = protocol.load_module('temperature module gen2', 'D1')
    temp_adapter = temp_block.load_adapter('opentrons_96_deep_well_adapter')
    shaker = protocol.load_module('heaterShakerModuleV1','C1')
    shake_adapter = shaker.load_adapter('opentrons_96_deep_well_adapter')
    mag_block = protocol.load_module('magneticBlockV1', 'D2')
    TRASH = protocol.load_waste_chute()
    EMPTYDECKSLOT      = '' # no empty deck slots
    EMPTYOFFDECKSLOT = 'A4' # for swapping labware
    TIP50_DECKSLOT = '' # none to start
    TIP50_OFFDECKSLOT = ''
    TIP1000_DECKSLOT = 'C3'
    TIP1000_OFFDECKSLOT = ''
    # endregion
    # region ================================ Helper Functions ================================
    global COLUMN_1_LIST
    COLUMN_1_LIST = ['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12']
    COUNTERS = {'WASTEVOL': 0}

    def swap_labware(labware1,labware2,EMPTYSLOT = '', gripper = True):
        labware1_pos = labware1.parent
        labware2_pos = labware2.parent
        if EMPTYSLOT:
            protocol.move_labware(labware1,EMPTYSLOT, use_gripper = gripper)
            protocol.move_labware(labware2,labware1_pos, use_gripper = gripper)
            protocol.move_labware(labware1,labware2_pos, use_gripper = gripper)
        else:
            protocol.move_labware(labware1,protocol_api.OFF_DECK)
            protocol.move_labware(labware2,labware1_pos, use_gripper = False)
            protocol.move_labware(labware1,labware2_pos, use_gripper = False)
    # endregion
    # region =============================== Tip Management ===========================
    ACTIVE_TIPRACKS = {"tip50": None, "tip1000": None}
    BACKUP_TIPRACKS = {"tip50": [], "tip1000": []}
    TIPS_USED = {"tip50": 0, "tip1000": 0} # tracks how many tips have been used in the active box

    def load_active_tips(active_tipracks, tips_used, tip_type, tip_rack):
        active_tipracks[tip_type] = tip_rack
        tips_used[tip_type] = 0
    def load_backup_tips(backup_tipracks, tip_type, tip_rack):
        backup_tipracks[tip_type].append(tip_rack)
    
    load_active_tips(ACTIVE_TIPRACKS, TIPS_USED, 'tip1000', protocol.load_labware(TIP1000_APINAME,TIP1000_DECKSLOT,TIP1000_APINAME))
    #load_active_tips(ACTIVE_TIPRACKS, TIPS_USED, 'tip50', protocol.load_labware(TIP50_APINAME,TIP50_OFFDECKSLOT,TIP50_APINAME))
    # load the backup slot
    #load_backup_tips(BACKUP_TIPRACKS,'tip1000',protocol.load_labware(TIP1000_APINAME,TIP1000_OFFDECKSLOT,TIP1000_APINAME))
        
    def load_new_tips(tip_type, tip_apiname, active_tipracks, backup_tipracks, tips_used, EMPTYOFFDECKSLOT):
        """
        Swap an empty tip box out for a new one from off deck.
        """
        protocol.comment("---> Swapping empty tips")
        old_box = active_tipracks[tip_type]
        old_box_location = old_box.parent
        if backup_tipracks[tip_type]: # if there are backup tips on deck
            new_box = backup_tipracks[tip_type].pop() # get the next box from the backup
            if EMPTYOFFDECKSLOT:
                swap_labware(old_box,new_box,EMPTYOFFDECKSLOT,gripper = True)
            else:
                swap_labware(old_box,new_box, gripper = False)
        else:
            new_box = protocol.load_labware(tip_apiname,protocol_api.OFF_DECK,tip_apiname)
            protocol.move_labware(old_box,protocol_api.OFF_DECK)
            protocol.move_labware(new_box,old_box_location)
        load_active_tips(active_tipracks, tips_used, tip_type, new_box)

    def get_next_tip(pipette, tip_type, tip_apiname, active_tipracks, backup_tipracks, tips_used):
        """
        Get the next tip from the tip rack.
        If the tips are empty, calls load_new_tips
        """
        if tips_used[tip_type] == 12:
            load_new_tips(tip_type, tip_apiname, active_tipracks, backup_tipracks, tips_used, EMPTYOFFDECKSLOT)
        tip_col = f"A{tips_used[tip_type] + 1}"
        pipette.pick_up_tip(active_tipracks[tip_type][tip_col])
        tips_used[tip_type] += 1

    # endregion
    # =============================== Define Labware ===============================
    #RES96_TYPE        = "vwr_96_wellplate_2500ul"
    RES96_TYPE       = "nest_96_wellplate_2ml_deep_custom"
    
    SamplePlate        = shake_adapter.load_labware(RES96_TYPE, 'Sample Plate')
    ReagentPlate        = protocol.load_labware(RES96_TYPE, 'C2', 'reagent reservoir')
    ElutionPlate    = temp_adapter.load_labware(RES96_TYPE,'Elution Plate')

    def transfer_tracktips(pipette, volume, source, dest, tip_type, tip_apiname, active_tiplist, backup_tiplist, tips_used, **kwargs):
        """
        Transfer liquid with tracking of tips.
        """
        # check if source is not a list
        if isinstance(source,list):
            if isinstance(source[0],list):
                raise ValueError("Source is more than one column")
        get_next_tip(pipette, tip_type, tip_apiname, active_tiplist, backup_tiplist, tips_used)
        pipette.transfer(volume, source, dest, new_tip = "never", **kwargs)
        pipette.drop_tip()

    def quench(Reagent,Sample,Elute,reagent_column,reagent_volume,sample_column,sample_volume):
        shaker.deactivate_shaker()
        transfer_tracktips(p1000, reagent_volume, Reagent[reagent_column], Sample[sample_column],
                            'tip1000', TIP1000_APINAME, ACTIVE_TIPRACKS, BACKUP_TIPRACKS, TIPS_USED)
        transfer_tracktips(p1000, sample_volume, Sample[sample_column], Elute[sample_column],
                            'tip1000', TIP1000_APINAME, ACTIVE_TIPRACKS, BACKUP_TIPRACKS, TIPS_USED,
                            mix_after=(5, sample_volume*0.75))
        shaker.set_and_wait_for_shake_speed(RPM)
        
    
    
    #region ================================ Protocol Steps ================================
    # SamplePlate starts on the shaker
    temp_block.set_temperature(4 if not DRYRUN else 25)
    if USETEMP:
        shaker.set_and_wait_for_temperature(TEMP)
    shaker.open_labware_latch()
    protocol.pause('Resume at t=0')
    shaker.close_labware_latch()
    shaker.set_and_wait_for_shake_speed(RPM)
    protocol.delay(seconds = 300 if not DRYRUN else 0.01)

    # t = 5
    quench(ReagentPlate,SamplePlate,ElutionPlate,'A1',100,'A1',1000)
    protocol.delay(seconds = 290 if not DRYRUN else 0.01)

    # t = 10
    quench(ReagentPlate,SamplePlate,ElutionPlate,'A1',100,'A2',1000)
    protocol.delay(seconds = 290 if not DRYRUN else 0.01)

    # t = 15
    quench(ReagentPlate,SamplePlate,ElutionPlate,'A1',100,'A3',1000)
    protocol.delay(seconds = 290 if not DRYRUN else 0.01)

    # t = 20
    quench(ReagentPlate,SamplePlate,ElutionPlate,'A1',100,'A4',1000)
    protocol.delay(seconds = 290 if not DRYRUN else 0.01)

    # t = 25
    quench(ReagentPlate,SamplePlate,ElutionPlate,'A1',100,'A5',1000)
    protocol.delay(seconds = 290 if not DRYRUN else 0.01)

    # t = 30
    quench(ReagentPlate,SamplePlate,ElutionPlate,'A1',100,'A6',1000)
    
    shaker.deactivate_shaker()
    shaker.open_labware_latch()
