# Overview of code:
'''
1. Take Parameters for dry run, number of columns, input volume, wash volumes, whether to include Dnase treatment, Dnase volume, and elution volume
3. will have two 96-well plates for ethanol washes. One for first wash, one for second wash
4. functions:
    wash(plate, columns, empty_space, wash_vol, wash_reservoir, pipette, tip_rack)

To simulate:
pixi run opentrons_simulate -L ./labware protocols/BombBio_Trizol_RNA_reuse_noShake_v3.py
'''

from opentrons import protocol_api
from opentrons import types
import math, re

metadata = {'protocolName': 'BombBio-trizol-rna-reuse-noShake-v3','author': 'Jared Bard','source': 'Protocol Library',}
requirements = {"robotType": "Flex","apiLevel": "2.21",}

def add_parameters(parameters):
    # ======================== RUNTIME PARAMETERS ========================
    parameters.add_bool(
        display_name="Dry Run",
        variable_name="DRYRUN",
        default=False,
        description="Whether to perform a dry run or not.")
    parameters.add_int(
        display_name="Sample Column count",
        variable_name="N_SAMPLECOLS",
        default=1,minimum=1,maximum=12,
        description="How many sample columns to process (start with 1).")
    parameters.add_float(
        display_name="Input Volume (µl)",
        variable_name="INPUTVOLUME",
        default=420,minimum=220,maximum=900,
        description="Input volume of sample + binding buffer + beads")
    parameters.add_float(
        display_name="EtOH Wash 1 Volume (µL)",
        variable_name="WASH1VOL",
        default=900,minimum=100,maximum=900,
        description="Volume for first wash")
    parameters.add_float(
        display_name="EtOH Wash 2 Volume (µL)",
        variable_name="WASH2VOL",
        default=500,minimum=100,maximum=900,
        description="Volume for subsequence EtOH wash")
    parameters.add_bool(
        display_name="Dnase Treat",
        variable_name="DNASETREAT",
        default=True,
        description="Include the Dnase Treatment")
    parameters.add_float(
        display_name="Dnase volume",
        variable_name="DNASEVOL",
        default=100,minimum=50,maximum=150,
        description="Volume for Dnase Treatment")
    parameters.add_float(
        display_name="Re-binding buffer volume",
        variable_name="REBINDVOL",
        default=500,minimum=100,maximum=1000,
        description="Re-binding buffer volume")
    parameters.add_float(
        display_name="Elution volume",
        variable_name="ELUTEVOL",
        default=50,minimum=20,maximum=500,
        description="Volume for final elution")
    

def run(protocol: protocol_api.ProtocolContext):
    # region ======================== DOWNLOADED PARAMETERS ========================
    DRYRUN              = protocol.params.DRYRUN
    N_SAMPLECOLS             = protocol.params.N_SAMPLECOLS
    INPUTVOLUME         = protocol.params.INPUTVOLUME
    WASH1VOL           = protocol.params.WASH1VOL
    WASH2VOL           = protocol.params.WASH2VOL
    DNASETREAT           = protocol.params.DNASETREAT
    DNASEVOL           = protocol.params.DNASEVOL
    REBINDVOL           = protocol.params.REBINDVOL
    ELUTEVOL           = protocol.params.ELUTEVOL

    # =================================================================================================
    # ====================================== ADVANCED PARAMETERS ======================================
    # =================================================================================================
    #---------------------------
    DNASETIME          = 15
    SHAKE_RPM           = 1300
    BINDTIME           = 10
    SETTLETIME          = 2
    DRYTIME          = 5
    EMPTYDECKSLOT      = 'C2'

    TIP_TRASH           = True      # True = Used tips go in Trash, False = Used tips go back into rack
    DEACTIVATE_TEMP     = True      # Whether or not to deactivate the heating and cooling modules after a run
    CUSTOM_OFFSETS      = False     # Manually enter offset position settings
    TIP50_apiname = 'opentrons_flex_96_filtertiprack_50ul'
    TIP1000_APINAME = 'opentrons_flex_96_filtertiprack_1000ul'
    TIP50_DECKSLOT = protocol_api.OFF_DECK
    TIP1000_DECKSLOT = 'B3'
    # endregion
    # region ============================ CUSTOM OFFSETS ===========================
    # These are Custom Offsets which are a PER INSTRUMENT Setting, to account for slight adjustments of the gripper calibration or labware.
    if CUSTOM_OFFSETS == True:
        PCRPlate_Z_offset = 0
        Deepwell_Z_offset = 1
        # HEATERSHAKER OFFSETS
        hs_drop_offset={'x':0,'y':-2,'z':0}
        hs_pick_up_offset={'x':0,'y':-2,'z':0}
        # MAG BLOCK OFFSETS
        mb_drop_offset={'x':0,'y':0.,'z':0.5}
        mb_pick_up_offset={'x':0,'y':-2,'z':0}
        # THERMOCYCLER OFFSETS
        tc_drop_offset={'x':0,'y':0,'z':0}
        tc_pick_up_offset={'x':0,'y':0,'z':0}
        # DECK OFFSETS
        deck_drop_offset={'x':0,'y':0,'z':0}
        deck_pick_up_offset={'x':0,'y':0,'z':0}
    else:
        PCRPlate_Z_offset = 0
        Deepwell_Z_offset = 0
        # HEATERSHAKER OFFSETS
        hs_drop_offset={'x':0,'y':0,'z':0}
        hs_pick_up_offset={'x':0,'y':0,'z':0}
        # MAG BLOCK OFFSETS
        mb_drop_offset={'x':0,'y':0.,'z':0}
        mb_pick_up_offset={'x':0,'y':0,'z':0}
        # THERMOCYCLER OFFSETS
        tc_drop_offset={'x':0,'y':0,'z':0}
        tc_pick_up_offset={'x':0,'y':0,'z':0}
        # DECK OFFSETS
        deck_drop_offset={'x':0,'y':0,'z':0}
        deck_pick_up_offset={'x':0,'y':0,'z':0}
    # endregion
    # region =============================== PIPETTE ===============================
    p50 = protocol.load_instrument("flex_8channel_50", "left")
    p1000 = protocol.load_instrument('flex_8channel_1000', 'right')
    p1000_flow_rate_aspirate_default = 200
    p1000_flow_rate_dispense_default = 200
    p1000_flow_rate_blow_out_default = 400
    p50_flow_rate_aspirate_default = 50
    p50_flow_rate_dispense_default = 50
    p50_flow_rate_blow_out_default = 100
    # endregion

    # region =============================== Deck Setup ===============================

    temp_block = protocol.load_module('temperature module gen2', 'C1')
    temp_adapter = temp_block.load_adapter('opentrons_96_well_aluminum_block')
    mag_block = protocol.load_module('magneticBlockV1', 'D1')
    TRASH = protocol.load_waste_chute()

    ACTIVE_TIPRACKS = {"tip50": {'deckslot': None, 'rack': None}, "tip1000": {'deckslot': TIP1000_DECKSLOT, 'rack': None}}
    BACKUP_TIPRACKS = {"tip50": [], "tip1000": []}
    ACTIVE_TIPRACKS['tip1000']['rack'] = protocol.load_labware(TIP1000_APINAME,TIP1000_DECKSLOT,"1000uL Filtered tips (1xUse)")
    # load the backup slot
    BACKUP_TIPRACKS['tip1000'].append(protocol.load_labware(TIP1000_APINAME,'D4',"1000uL Filtered tips (1xUse)"))
    TIPS_USED = {"tip50": 0, "tip1000": 0} # tracks how many tips have been used in the active box
    tip1000_reuse = protocol.load_labware(TIP1000_APINAME, 'A3')
        
    # =============================== Define Labware ===============================
    #RES96_TYPE        = "vwr_96_wellplate_2500ul"
    RES96_TYPE       = "nest_96_wellplate_2ml_deep"
    RES1_TYPE         = "nest_1_reservoir_290ml"
    PLATE_TYPE      =   "nest_96_wellplate_2ml_deep"
    REAGENT_PLATE_TYPE = "nest_96_wellplate_2ml_deep"
    ELUTION_PLATE_TYPE = "opentrons_96_wellplate_200ul_pcr_full_skirt"

    SamplePlate        = protocol.load_labware(RES96_TYPE, protocol_api.OFF_DECK, 'Sample Plate')
    Wash1Res        = protocol.load_labware(RES96_TYPE, 'A2', 'Wash1 reservoir EtOH')
    Wash2Res       = protocol.load_labware(RES96_TYPE, 'B2', 'Wash2 reservoir EtOH')
    WasteRes       = protocol.load_labware(RES1_TYPE, 'C3', 'Waste Reservoir')
    ReagentPlate        = protocol.load_labware(REAGENT_PLATE_TYPE, 'D2', 'reagent reservoir 3')
    ElutionPlate    = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt',protocol_api.OFF_DECK,'Elution Plate')
    
    # ======== DEFINING LIQUIDS =======
    SampleLiq = protocol.define_liquid(name="Sample", description="Sample", display_color="#E69F00")  # Orange
    EtOHLiq = protocol.define_liquid(name="EtOH", description="90% Ethanol", display_color="#56B4E9")  # Sky Blue
    DnaseLiq = protocol.define_liquid(name="DNase_Sol", description="DnaseI in buffer", display_color="#009E73")  # Bluish Green
    BindingBufferLiq = protocol.define_liquid(name="BindingBuffer", description="Binding Buffer", display_color="#F0E442")  # Yellow
    ElutionBufferLiq = protocol.define_liquid(name="ElutionBuffer", description="Elution buffer", display_color="#D55E00")  # Vermillion

    # ======== ESTIMATING LIQUIDS =======
    Sample_Volume = INPUTVOLUME
    WASH1_Vol_Per_Well = 2200
    WASH2_Vol_Per_Well = 2200
    DNASE_Vol_Per_Well = DNASEVOL * N_SAMPLECOLS * 1.1
    REBIND_Vol_Per_Well = 2200
    ELUTE_Vol_Per_Well = ELUTEVOL * N_SAMPLECOLS * 1.1

    # =============================== Load Liquids ===============================

    def load_column_liquid(labware, column_index, liq, vol):
        row_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        well_list = [f"{row}{column_index}" for row in row_labels]
        for well in well_list:
            protocol.comment(f"Loading liquid into {labware.wells_by_name()}")
            labware.wells_by_name()[well].load_liquid(liquid=liq, volume=vol)
    
    sample_column_keys = [str(i) for i in range(1, N_SAMPLECOLS+1)]
    Sample_cols = {key: SamplePlate.columns_by_name()[key] for key in sample_column_keys}
    Wash1EtOH_cols = {key: Wash1Res.columns_by_name()[key] for key in sample_column_keys}
    Wash2EtOH_cols = {key: Wash2Res.columns_by_name()[key] for key in sample_column_keys}

    DnaseBuffer_cols = ReagentPlate.columns_by_name()['1']  # first column

    if N_SAMPLECOLS <= 4:
        binding_buffer_column = [3]
    elif N_SAMPLECOLS > 4 and N_SAMPLECOLS <= 8:
        binding_buffer_column = [3, 4]
    elif N_SAMPLECOLS > 8 and N_SAMPLECOLS <= 12:
        binding_buffer_column = [3, 4, 5]
    else:
        raise ValueError("Unsupported number of samples")
    
    BindingBuffer_cols = {key: ReagentPlate.columns_by_name()[key] for key in [str(i) for i in binding_buffer_column]}
    ElutionBuffer_cols = ReagentPlate.columns_by_name()['12']  # last column

    # Grouped load_column_liquid commands
    for i in range(1, N_SAMPLECOLS+1):
        load_column_liquid(SamplePlate, i, SampleLiq, Sample_Volume)
        load_column_liquid(Wash1Res, i, EtOHLiq, WASH1_Vol_Per_Well)
        load_column_liquid(Wash2Res, i, EtOHLiq, WASH2_Vol_Per_Well)

    load_column_liquid(ReagentPlate, 1, DnaseLiq, DNASE_Vol_Per_Well)

    for i in binding_buffer_column:
        load_column_liquid(ReagentPlate, i, BindingBufferLiq, REBIND_Vol_Per_Well)

    load_column_liquid(ReagentPlate, 1, ElutionBufferLiq, ELUTE_Vol_Per_Well)
    # endregion
    # region ================================ Helper Functions ================================
    global COLUMN_1_LIST
    COLUMN_1_LIST = ['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12']
    SAMPLECOLS = [COLUMN_1_LIST[i] for i in range(N_SAMPLECOLS)]
    COUNTERS = {'WASTEVOL': 0}
    # endregion
    
    def swap_tips(tip_type, tip_apiname, active_tipracks, backup_tipracks, new_slot = None):
        """
        Swaps an empty tip box out for a new one from off deck.
        First define a tip box by iterating 1 from the offdeck_count
        Then swap the tip box in from off_deck
        Returns the new tipbox_count
        """
        old_box = active_tipracks[tip_type]['rack']
        # Extract the digits at the end of the variable name using regex
        if new_slot:
            active_tipracks[tip_type]['deckslot'] = new_slot
        if not old_box: # if there is no rack on deck
            new_box = protocol.load_labware(tip_apiname, protocol_api.OFF_DECK)
            protocol.move_labware(labware = new_box, new_location = active_tipracks[tip_type]['deckslot'])
        else:
            if backup_tipracks[tip_type]: # if there are backup tips on deck
                protocol.comment("---> Remove old tips")
                protocol.move_labware(labware = old_box, new_location = 'A4', use_gripper=True) # drop off the deck
                protocol.comment("---> Grabbing backup tips from on-deck")
                new_box = backup_tipracks[tip_type].pop() # get the next box from the backup
                protocol.move_labware(labware = new_box, new_location = active_tipracks[tip_type]['deckslot'],
                              use_gripper=True)
            else:
                protocol.comment("---> Remove old tips")
                new_box = protocol.load_labware(tip_apiname, protocol_api.OFF_DECK)
                protocol.move_labware(labware = old_box, new_location = protocol_api.OFF_DECK) # drop off the deck
                protocol.comment("---> Swap backup tips from off-deck")
                protocol.move_labware(labware = new_box, new_location = active_tipracks[tip_type]['deckslot'])
        active_tipracks[tip_type]['rack'] = new_box

    def get_next_tip(pipette, tip_type, tip_apiname, active_tipracks, backup_tipracks, tips_used):
        """
        Get the next tip from the tip rack.
        If the tips are empty, calls swap_tips
        """
        if tips_used[tip_type] == 12:
            swap_tips(tip_type, tip_apiname, active_tipracks, backup_tipracks)
            tips_used[tip_type] = 0
        tip_col = f"A{tips_used[tip_type] + 1}"
        pipette.pick_up_tip(active_tipracks[tip_type]['rack'][tip_col])
        tips_used[tip_type] += 1

    def transfer_reuse(pipette, volume, source, dest, tips, **kwargs):
        """
        Wrapper function for transferring liquid with reusable tips.
        
        Args:
            pipette: The pipette object (e.g., an instance of InstrumentContext).
            volume: The volume to transfer (float or Sequence[float]).
            source: The source location (AdvancedLiquidHandling).
            dest: The destination location (AdvancedLiquidHandling).
            trash: Whether to discard the tip after use. Defaults to True.
            **kwargs: Additional arguments to pass to the transfer method.
        
        Returns:
            The result of the transfer method.
        """
        pipette.pick_up_tip(tips)
        pipette.transfer(volume, source, dest, new_tip = "never", **kwargs)
        pipette.drop_tip(tips)

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

    def removeSup(PLATE, COL, SUPVOL, counter_dict, Deepwell_Z_offset = 0):
        """
        Remove supernatant from the plate using the p1000
        no tip tracking
        """
        if SUPVOL < 100:
            SUPVOL = 100
        p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default*0.5
        p1000.flow_rate.dispense = p1000_flow_rate_dispense_default*0.5
        p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default*0.5
        p1000.move_to(PLATE[COL].bottom(z=Deepwell_Z_offset+2))
        p1000.aspirate(SUPVOL-100)
        protocol.delay(minutes=0.1)
        p1000.move_to(PLATE[COL].bottom(z=Deepwell_Z_offset+1))
        p1000.aspirate(100)
        p1000.default_speed = 200
        p1000.move_to(PLATE[COL].top(z=2))
        #======L Waste Volume Check======
        counter_dict['WASTEVOL'] += (SUPVOL*8)
        protocol.comment('--->Adding '+str((SUPVOL*8)/1000)+'mL to waste (total is '+str(counter_dict['WASTEVOL']/1000)+'mL)')
        if counter_dict['WASTEVOL'] >150000:
            protocol.pause('Please empty the waste')
        p1000.dispense(SUPVOL, WasteRes['A1'].top(z=0))
        protocol.delay(minutes=0.1)
        p1000.blow_out()
        p1000.default_speed = 400
        p1000.move_to(WasteRes['A1'].top(z=-5))
        p1000.move_to(WasteRes['A1'].top(z=0))
        p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default
        p1000.flow_rate.dispense = p1000_flow_rate_dispense_default
        p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default
    
    def wash_plate(SAMPLEPLATE, SAMPLEPOS, COLS, WASHVOL, WASHRES, SETTLETIME, tip_rack, counter_dict, Deepwell_Z_offset):
        """
        Wash the plate with EtOH.
        Should change this to start with an empty plate on the magblock
        """

        protocol.comment('--> Adding Wash')
        for i, X in enumerate(COLS):
            p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default*0.5
            p1000.flow_rate.dispense = p1000_flow_rate_dispense_default
            p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default
            transfer_reuse(p1000, WASHVOL, WASHRES[X].bottom(Deepwell_Z_offset), SAMPLEPLATE[X], tip_rack[X], mix_after=(3, WASHVOL*0.75),
                           air_gap = 20, blow_out=True, blowout_location="trash")

        p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default
        p1000.flow_rate.dispense = p1000_flow_rate_dispense_default
        p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default
        protocol.comment('--> Moving Sample Plate to MagBlock')
        protocol.move_labware(labware=SAMPLEPLATE,new_location=mag_block,use_gripper=True)

        protocol.comment('--> Wait for beads to settle')
        if DRYRUN == False:
            protocol.delay(minutes=SETTLETIME)
            
        protocol.comment('--> Remove Supernatant')
        for i, X in enumerate(COLS):
            p1000.pick_up_tip(tip_rack[X])
            removeSup(SAMPLEPLATE, X, WASHVOL+25, counter_dict, Deepwell_Z_offset)
            p1000.drop_tip(tip_rack[X])
        
        protocol.comment('--> Moving plate off magnet')
        protocol.move_labware(labware = SAMPLEPLATE, new_location = SAMPLEPOS, use_gripper = True)

    #region ================================ Protocol Steps ================================
    # start with the sample plate on the magblock
    protocol.move_labware(labware=SamplePlate,new_location=mag_block)

    protocol.comment('--> Wait for beads to settle')
    if DRYRUN == False:
        protocol.delay(minutes=SETTLETIME)
        
    protocol.comment('--> Remove Sample')
    for i, X in enumerate(SAMPLECOLS):
        get_next_tip(p1000, 'tip1000', TIP1000_APINAME, ACTIVE_TIPRACKS, BACKUP_TIPRACKS, TIPS_USED)
        removeSup(SamplePlate, X, INPUTVOLUME, COUNTERS, Deepwell_Z_offset)
        p1000.drop_tip()
    
    protocol.comment('--> Moving plate off magnet')
    protocol.move_labware(labware=SamplePlate,new_location=protocol_api.OFF_DECK)
    protocol.pause('Spin plate at 500g for 1 minute to collect residual trizol. meanwhile, empty the trizol waste into the waste reservoir and rinse the waste container')
    protocol.move_labware(labware=SamplePlate,new_location=EMPTYDECKSLOT)

    protocol.comment('--> Wash1')
    wash_plate(SamplePlate, EMPTYDECKSLOT, SAMPLECOLS, WASH1VOL, Wash1Res, 2, tip1000_reuse, COUNTERS, 20)
    protocol.comment('--> Wash2')
    wash_plate(SamplePlate, EMPTYDECKSLOT, SAMPLECOLS, WASH2VOL, Wash2Res, 2, tip1000_reuse, COUNTERS, 20)
    protocol.comment('--> Wash3')
    wash_plate(SamplePlate, EMPTYDECKSLOT, SAMPLECOLS, WASH2VOL, Wash2Res, 2, tip1000_reuse, COUNTERS, Deepwell_Z_offset)


    protocol.comment('--> Dry Plate')
    if DRYRUN == False:
        temp_block.set_temperature(50)
        temp_block.await_temperature(50)
    else:
        temp_block.set_temperature(50)
        temp_block.await_temperature(50)
    protocol.move_labware(labware=SamplePlate,new_location=temp_adapter, use_gripper=True)
    if DRYRUN == False:
        protocol.delay(minutes=DRYTIME)
    if DRYRUN == False:
        temp_block.set_temperature(37)
        temp_block.await_temperature(37)
    else:
        temp_block.set_temperature(37)
    
    protocol.comment('--> Add DNaseI')
    for X in SAMPLECOLS:
        transfer_tracktips(p1000, DNASEVOL, DnaseBuffer_cols, SamplePlate[X],
        'tip1000', TIP1000_APINAME, ACTIVE_TIPRACKS, BACKUP_TIPRACKS, TIPS_USED,
        mix_after=(10, DNASEVOL*0.75))
    
    protocol.comment('--> Incubate')
    if DRYRUN == False:
        protocol.delay(minutes=DNASETIME)
    
    protocol.comment('--> Rebind')
    # use the same tip to distribute binding buffer to each column
    get_next_tip(p1000, 'tip1000', TIP1000_APINAME, ACTIVE_TIPRACKS, BACKUP_TIPRACKS, TIPS_USED)
    p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default*0.5
    p1000.flow_rate.dispense = p1000_flow_rate_dispense_default*0.5
    p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default*0.5
    for i, X in enumerate(BindingBuffer_cols.keys()):
        source_well = ReagentPlate[f"A{X}"]
        # only 4 transfers per column of binding buffer, so switch to next column after 4
        # multidispense the binding buffer because we will shake to mix afterwards anyways
        for Y in SAMPLECOLS[i*4:i*4+3]:
            p1000.aspirate(REBINDVOL, source_well.bottom(z=1))
            p1000.air_gap(20)
            p1000.move_to(source_well.top(z=0))
            p1000.move_to(source_well.top(z=-5))
            p1000.move_to(source_well.top(z=0))
            #=====Reservoir Tip Touch========
            p1000.default_speed = 100
            p1000.move_to(source_well.top().move(types.Point(x=4,z=-3)))
            p1000.move_to(source_well.top().move(types.Point(x=-4,z=-3)))
            p1000.default_speed = 400
            #================================ 
            p1000.move_to(SamplePlate[Y].top(z=7))
            p1000.dispense(REBINDVOL+20)
            protocol.delay(minutes=0.1)
            p1000.move_to(SamplePlate[Y].top(z=5))
            p1000.move_to(SamplePlate[Y].top(z=2))
            p1000.move_to(SamplePlate[Y].top(z=5))
            p1000.air_gap(20) # to prevent leaking while moving
    p1000.drop_tip()
    p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default
    p1000.flow_rate.dispense = p1000_flow_rate_dispense_default
    p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default
    
    protocol.comment('--> Shake')
    protocol.move_labware(labware=SamplePlate,new_location=protocol_api.OFF_DECK)
    protocol.pause('Shake for 10 minutes at 1300rpm')
    protocol.move_labware(labware=SamplePlate,new_location=mag_block)
    protocol.comment('--> Wait for beads to settle')
    if DRYRUN == False:
        protocol.delay(minutes=SETTLETIME)
        
    protocol.comment('--> Remove Sample')
    for i, X in enumerate(SAMPLECOLS):
        get_next_tip(p1000, 'tip1000', TIP1000_APINAME, ACTIVE_TIPRACKS, BACKUP_TIPRACKS, TIPS_USED)
        removeSup(SamplePlate, X, REBINDVOL, COUNTERS, Deepwell_Z_offset)
        p1000.drop_tip()
    
    protocol.comment('--> Moving plate off magnet')
    protocol.move_labware(labware=SamplePlate,new_location=protocol_api.OFF_DECK)
    protocol.pause('Spin plate at 500g for 1 minute to collect residual trizol. meanwhile, empty the trizol waste into the waste reservoir and rinse the waste container')
    protocol.move_labware(labware=SamplePlate,new_location=EMPTYDECKSLOT)
    
    protocol.comment('--> Wash1')
    wash_plate(SamplePlate, EMPTYDECKSLOT, SAMPLECOLS, WASH1VOL, Wash1Res, 2, tip1000_reuse, COUNTERS, Deepwell_Z_offset)
    protocol.comment('--> Wash2')
    wash_plate(SamplePlate, EMPTYDECKSLOT, SAMPLECOLS, WASH2VOL, Wash2Res, 2, tip1000_reuse, COUNTERS, Deepwell_Z_offset)
    protocol.comment('--> Wash3')
    wash_plate(SamplePlate, EMPTYDECKSLOT, SAMPLECOLS, WASH2VOL, Wash2Res, 2, tip1000_reuse, COUNTERS, Deepwell_Z_offset)
    
    protocol.comment('--> Disposing of tips')
    for X in SAMPLECOLS:
        p1000.pick_up_tip(tip1000_reuse[X])
        p1000.drop_tip()

    protocol.comment('--> Dry Plate')
    if DRYRUN == False:
        temp_block.set_temperature(50)
        temp_block.await_temperature(50)
    else:
        temp_block.set_temperature(50)

    protocol.move_labware(labware=SamplePlate,new_location=temp_adapter, use_gripper=True)
    if DRYRUN == False:
        protocol.delay(minutes=DRYTIME)
    protocol.move_labware(labware=SamplePlate,new_location=EMPTYDECKSLOT, use_gripper=True)

    protocol.comment('--> Elute')
    # load the p50 tips
    protocol.move_labware(labware = tip1000_reuse, new_location=protocol_api.OFF_DECK)
    swap_tips('tip50', TIP50_apiname, ACTIVE_TIPRACKS, BACKUP_TIPRACKS, new_slot = 'A3')
    protocol.move_labware(labware=Wash2Res,new_location=protocol_api.OFF_DECK)
    protocol.move_labware(labware=ElutionPlate,new_location='B2')
    for i, X in enumerate(SAMPLECOLS):
        transfer_tracktips(p50, ELUTEVOL, ElutionBuffer_cols, SamplePlate[X],
        'tip50', TIP50_apiname, ACTIVE_TIPRACKS, BACKUP_TIPRACKS, TIPS_USED,
        mix_after=(10, ELUTEVOL*0.75))
    if DRYRUN == False:
        protocol.delay(minutes=5)
    protocol.move_labware(labware = SamplePlate, new_location = mag_block, use_gripper = True)
    protocol.comment('--> Wait for beads to settle')
    if DRYRUN == False:
        protocol.delay(minutes=SETTLETIME)
    protocol.comment('--> Elute')
    for i, X in enumerate(SAMPLECOLS):
        transfer_tracktips(p50, ELUTEVOL, SamplePlate[X], ElutionPlate[X],
                           'tip50', TIP50_apiname, ACTIVE_TIPRACKS, BACKUP_TIPRACKS, TIPS_USED)