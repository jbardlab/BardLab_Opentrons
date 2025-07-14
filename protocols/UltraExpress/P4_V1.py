from opentrons import protocol_api
from opentrons import types

# from opentrons.protocol_api import SINGLE, ALL
import math
import numpy as np

metadata = {
    "protocolName": "NEBNext UltraExpress® RNA Library Prep Kit_NEB #E3330S/L_Part4_V1",
    "author": "Opentrons",
    "description": "Part 4:Phased Bead Cleanup of PCR Reaction",
}

requirements = {"robotType": "Flex", "apiLevel": "2.20"}


def add_parameters(parameters):
    parameters.add_bool(
        variable_name="DRY_WATER_RUN",
        display_name="Dry or Sample Run",
        description="Do you want to perform a dry run?",
        default=False,
    )
    parameters.add_int(
        variable_name="num_samples",
        display_name="Number of Samples",
        description="The number of samples to process",
        default=94,
        minimum=1,
        maximum=96,
        unit="x Samples",
    )
    parameters.add_str(
        variable_name="p50m_mount",
        display_name="50 µL Pipette Mount",
        description="The mount to attach the 50 µL pipette to",
        choices=[
            {"display_name": "Left", "value": "left"},
            {"display_name": "Right", "value": "right"},
        ],
        default="left",
    )
    parameters.add_str(
        variable_name="p1000m_mount",
        display_name="1000 µL Pipette Mount",
        description="The mount to attach the 1000 µL pipette to",
        choices=[
            {"display_name": "Left", "value": "left"},
            {"display_name": "Right", "value": "right"},
        ],
        default="right",
    )
    parameters.add_float(
        variable_name="flow_rate_aspirate",
        display_name="Default Flow rate aspirate",
        description="Flow rate for the aspirate step",
        default=700,
        minimum=5,
        maximum=716,
        unit="µL/s",
    )
    parameters.add_float(
        variable_name="flow_rate_dispense",
        display_name="Default Flow rate dispense",
        description="Flow rate for the dispense step",
        default=700,
        minimum=5,
        maximum=716,
        unit="µL/s",
    )
    parameters.add_float(
        variable_name="sample_rate",
        display_name="Sample Rate",
        description="Sample mixing rate",
        default=0.6,
        minimum=0,
        maximum=1.5,
    )
    parameters.add_float(
        variable_name="beads_rate",
        display_name="Beads Rate",
        description="Beads mixing rate",
        default=0.2,
        minimum=0,
        maximum=1.5,
    )
    parameters.add_float(
        variable_name="buffer_rate",
        display_name="Buffer Rate",
        description="Buffer mixing rate",
        default=0.8,
        minimum=0,
        maximum=1.5,
    )
    parameters.add_float(
        variable_name="ethanol_rate",
        display_name="Ethanol Rate",
        description="Ethanol mixing rate",
        default=0.8,
        minimum=0,
        maximum=1.5,
    )
    parameters.add_float(
        variable_name="elution_rate",
        display_name="Elution Rate",
        description="Elution mixing rate",
        default=0.4,
        minimum=0,
        maximum=1.5,
    )


def run(protocol: protocol_api.ProtocolContext):
    # ------------------------------------------------------------------ #
    #                             Parameters                             #
    # ------------------------------------------------------------------ #

    # ------------------ Load the run time parameters ------------------ #
    DRY_WATER_RUN = protocol.params.DRY_WATER_RUN
    num_samples = protocol.params.num_samples
    p50m_mount = protocol.params.p50m_mount
    p1000m_mount = protocol.params.p1000m_mount
    flow_rate_aspirate = protocol.params.flow_rate_aspirate
    flow_rate_dispense = protocol.params.flow_rate_dispense
    sample_rate = protocol.params.sample_rate
    beads_rate = protocol.params.beads_rate
    buffer_rate = protocol.params.buffer_rate
    ethanol_rate = protocol.params.ethanol_rate
    elution_rate = protocol.params.elution_rate
    # ------------------------ volume parameters ----------------------- #
    sample_vol = 100
    res_beads_vol = 0.7 * sample_vol
    te_1_vol = 50
    recon_buffer_vol = 0.8 * te_1_vol
    ethanol_vol = 200
    te_2_in_vol = 23
    te_2_out_vol = 20

    # -------------------- magnetic block parameters ------------------- #
    beads_incubation_time = 2
    beads_incubation_time_long = 5
    beads_binding_time = 5
    air_dry_time = 5

    # ------------------------------------------------------------------ #
    #                           Error Handling                           #
    # ------------------------------------------------------------------ #

    # 1. Check if the number of samples is within the range
    if num_samples > 96 or num_samples < 1:
        raise ValueError("The number of samples should be between 1 and 96")
    # 2. Any of the volumes should be equal to or less than 200 µL
    if (
        sample_vol > 200
        or res_beads_vol > 200
        or te_1_vol > 200
        or recon_buffer_vol > 200
        or ethanol_vol > 200
        or te_2_in_vol > 200
        or te_2_out_vol > 200
    ):
        raise ValueError("The volume should be equal to or less than 200 µL")

    # ------------------------------------------------------------------ #
    #                             Calculation                            #
    # ------------------------------------------------------------------ #
    # 1. Calculate the column number
    cols = math.ceil(num_samples / 8)
    if num_samples % 8 == 0:
        cols_m = cols
        cols_s = None
    elif num_samples % 8 != 0:
        if cols == 1:
            cols_m = 0
            cols_s = num_samples
        else:
            cols_m = cols - 1
            cols_s = num_samples % 8
    # print(f"cols: {cols}, cols_m: {cols_m}, cols_s: {cols_s}")
    # cols_m and cols_s will be used in the dt25_beads preparation step

    # ------------------------------------------------------------------ #
    #                Tips handling calculation and loading               #
    # ------------------------------------------------------------------ #

    # tip usage calculation
    # Rule 1: Use 50 µL tips for elution steps.
    # Rule 2: Other steps use 200 µL tips for better beads mixing and removal.

    # --------------------- step by step analysis --------------------- #
    # ------------------- calculate 200 ul tips first ------------------ #
    # 1. steps 1A9.1 to 1A9.5: 70 µL beads in, 170 µL supernatant out
    # cols x column for beads addition,
    # cols x column for supernatant removal
    # move to the top of temperature module for the next step
    required_cols_200_1 = cols + cols

    # 2. steps 1A9.6 to 1A9.9: 50 µL TE in, 40 µL NEBNext Bead Reconstitution Buffer in, 90 uL supernatant out
    # cols x column for TE addition and mixing
    # cols x column for buffer addition and mixing
    # cols x column for supernatant removal
    # plate stay on the magnetic block for ethanol wash
    # Use 200 µL tips for better mixing and removal
    required_cols_200_2 = cols + cols + cols

    # 3. steps 1A9.10: 200 µL 80% Ethanol in, 200 µL 80% Ethanol out
    # 1 column for ethanol addition on the top of each well, 30 seconds incubation
    # 1 column for supernatant removal
    # no air dry for the 1st round of ethanol wash
    # plate stay on the magnetic block for the 2nd round of ethanol wash
    required_cols_200_3 = 1 + cols
    # required_cols_200_2 = cols + cols

    # 4. steps 1A9.11 to 1A9.12: 200 µL 80% Ethanol in, 200 µL 80% Ethanol out
    # 1 column for ethanol addition on the top of each well, 30 seconds incubation
    # 1 column for supernatant removal
    # cols x columns: supernatant removal
    # air dry for 5 minutes
    # move plate to the top of the temperature module for the next step
    required_cols_200_4 = 1 + cols
    # required_cols_200_3 = cols + cols

    # ------------------- calculate 50 ul tips first ------------------- #
    # 6. steps 1A5.8: 23 µL TE in, 20 µL TE out
    # cols x columns: Adding TE buffer, mixing
    # cols x columns: cDNA removal (quickly)
    # new sample_plate should be kept on the thermocycler for the next step
    required_cols_50_1 = cols + cols

    required_cols_50 = required_cols_50_1
    required_slots_num_50 = math.ceil(required_cols_50 / 12)

    slots_50_extra = []
    expansion_slots_200 = ["A4", "B4", "C4", "D4"]
    if required_slots_num_50 >= 2:
        slots_50 = ["C3", "D2"]
        regular_slot_list = ["A2", "B2", "B3", "A3"]
        for i in range(required_slots_num_50 - 2):
            slots_50_extra.append(protocol_api.OFF_DECK)
    else:
        required_slots_num_50 = 1
        slots_50 = ["C3"]
        regular_slot_list = ["A2", "B2", "B3", "A3", "D2"]

    # ---------------------- calculate 200 ul tips --------------------- #
    required_cols_200 = (
        required_cols_200_1
        + required_cols_200_2
        + required_cols_200_3
        + required_cols_200_4
    )
    required_slots_num_200 = math.ceil(required_cols_200 / 12)

    # ------------------- assign slots for 200ul tips ------------------ #
    difference = required_slots_num_200 - len(regular_slot_list)
    if difference <= 0:
        slots_200 = regular_slot_list[:required_slots_num_200]
        slots_200_extra = []
        slots_empty_expansion = []
    else:
        slots_200 = regular_slot_list
        slots_200_extra = expansion_slots_200[:difference]
        slots_empty_expansion = expansion_slots_200[difference:]

    print(f"Number of samples: {num_samples} x samples")
    print(f"50ul tips number: {required_slots_num_50} x slots")
    print(f"50ul tips columns: {int(cols*2)}")
    print(f"50ul tips slots: {slots_50}")
    print(f"50ul tips extra slots: {slots_50_extra}")
    print(f"200ul tips number: {required_slots_num_200} x slots")
    print(f"200ul tips columns: {int(required_cols_200)}")
    print(f"200ul tips slots: {slots_200}")
    print(f"200ul tips extra slots: {slots_200_extra}")
    print(f"Slots not used: {slots_empty_expansion}")

    # ---------------------------- load tips --------------------------- #
    tips50 = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_50ul", slot)
        for slot in slots_50
    ]
    tips50_extra = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_50ul", slot)
        for slot in slots_50_extra
    ]
    tips200 = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_200ul", slot)
        for slot in slots_200
    ]
    tips200_extra = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_200ul", slot)
        for slot in slots_200_extra
    ]

    # ------------------------------------------------------------------ #
    #                               Pipettes                             #
    # ------------------------------------------------------------------ #
    p50m = protocol.load_instrument(
        "flex_8channel_50",
        p50m_mount,
        tips50,
    )
    p1000m = protocol.load_instrument(
        "flex_8channel_1000",
        p1000m_mount,
        tips200,
    )
    p1000m.flow_rate.aspirate = flow_rate_aspirate
    p1000m.flow_rate.dispense = flow_rate_dispense
    # # ------------------------------------------------------------------ #
    # #                               Modules                              #
    # # ------------------------------------------------------------------ #
    thermocycler = protocol.load_module("thermocycler module gen2")

    mag = protocol.load_module("magneticBlockV1", "D1")
    temp_mod = protocol.load_module("temperature module gen2", "C1")
    temp_adapter = temp_mod.load_adapter("opentrons_96_well_aluminum_block")
    chute = protocol.load_waste_chute()
    # # ------------------------------------------------------------------ #
    # #                               labware                              #
    # # ------------------------------------------------------------------ #
    reagent_plate = protocol.load_labware(
        "nest_96_wellplate_2ml_deep", "C2", "Reagent Plate"
    )
    sample_plate = temp_adapter.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt", "Sample Plate"
    )
    final_plate = thermocycler.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt",
        label="Final Plate",
    )
    # thermocycler.set_block_temperature(temperature=4)
    thermocycler.close_lid()
    # waste_res = protocol.load_labware("nest_1_reservoir_195ml", "D2", "Waste Reservoir")
    transition_slot = temp_adapter
    # ------------------------------------------------------------------ #
    #                           Plate mapping                           #
    # ------------------------------------------------------------------ #

    # _cols represent the 1st well of each column that will be used in the washing and elution steps
    # _cols_m represent the 1st well of each column that will be used in the full column preparation
    # _cols_s represent the leftover wells that will be used with single tip pickup
    # _wells represent all the wells that in sample_plate, used in the liquid labeling
    # _list represent the list of all the wells(from reagent plate or rt reagent plate) that will be used in the liquid labelling

    # -------------------------- Sample plate -------------------------- #
    sample_plate_cols = sample_plate.rows()[0][:cols]
    sample_plate_cols_m = sample_plate.rows()[0][:cols_m]
    if cols_s:
        sample_plate_cols_s = sample_plate.columns()[cols_m][:cols_s]
    else:
        sample_plate_cols_s = None
    sample_plate_wells = sample_plate.wells()[:num_samples]

    # print(f"sample_plate_cols: {sample_plate_cols}")
    # print(f"sample_plate_cols_m: {sample_plate_cols_m}")
    # print(f"sample_plate_cols_s: {sample_plate_cols_s}")
    # -------------------------- Reagent Plate ------------------------- #
    res_beads = reagent_plate.rows()[0][0]  # A1
    res_beads_list = reagent_plate.columns()[0]
    te = reagent_plate.rows()[0][1]  # A2
    te_list = reagent_plate.columns()[1]
    recon_buffer = reagent_plate.rows()[0][2]  # A3
    recon_buffer_list = reagent_plate.columns()[2]
    ethanol_1 = reagent_plate.rows()[0][3:5]  # A4 and A5, RT, 1st round
    ethanol_1_list = reagent_plate.columns()[3] + reagent_plate.columns()[4]
    ethanol_2 = reagent_plate.rows()[0][5:7]  # A6 and A7, RT, 2nd round
    ethanol_2_list = reagent_plate.columns()[5] + reagent_plate.columns()[6]
    # ------------------------- waste reservoir ------------------------ #
    waste_well = reagent_plate.rows()[0][8:12]  # A9 to A12, Waste

    # summarize reagent plate map in a list
    reagent_plate_list = [
        "Column 1 (A1 of Reagent Plate): Resuspension Beads, 70 µL per sample",
        "Column 2 (A2 of Reagent Plate): TE Buffer, 50 µL first for Reconstitution and 23 µL second for Elution",
        "Column 3 (A3 of Reagent Plate): NEBNext Bead Reconstitution Buffer, 40 µL per sample",
        "Column 4-5 (A4-A5 of Reagent Plate): 80% Ethanol, 1st round",
        "Column 6-7 (A6-A7 of Reagent Plate): 80% Ethanol, 2nd round",
        "Column 9-12(A9-A12 of Reagent Plate): Waste",
    ]
    print(reagent_plate_list)
    # ----------------------- Final Plate Part 1 ----------------------- #
    final_plate_cols = final_plate.rows()[0][:cols]
    final_plate_wells = final_plate.wells()[:num_samples]

    # ------------------------------------------------------------------ #
    #                               liquids                              #
    # ------------------------------------------------------------------ #
    l_locations = [
        res_beads_list,
        te_list,
        recon_buffer_list,
        ethanol_1_list if cols > 8 else ethanol_1_list[:8],
        ethanol_2_list if cols > 8 else ethanol_2_list[:8],
    ]
    l_volumes = [
        res_beads_vol * cols * 1.2,
        (te_1_vol + te_2_in_vol) * cols * 1.2,
        recon_buffer_vol * cols * 1.2,
        ethanol_vol * 8 * 1.2 if cols > 8 else ethanol_vol * cols * 1.2,
        ethanol_vol * 8 * 1.2 if cols > 8 else ethanol_vol * cols * 1.2,
    ]
    liquids = [
        "Resuspension Beads",
        "TE Buffer",
        "NEBNext Bead Reconstitution Buffer",
        "80% Ethanol",
        "80% Ethanol",
    ]
    descriptions = [
        "Resuspension Beads, 70 µL per sample",
        "TE Buffer, 50 µL per sample + 23 µL per sample",
        "NEBNext Bead Reconstitution Buffer, 40 µL per sample",
        "80% Ethanol, 200 µL per sample, 1st round, up to 2 columns",
        "80% Ethanol, 200 µL per sample, 2nd round, up to 2 columns",
    ]

    colors_full = [
        "#FF0000",  # Red
        "#0000FF",  # Blue
        "#008000",  # Green
        "#FFFF00",  # Yellow
        "#FFC0CB",  # Pink
        "#800080",  # Purple
        "#FFA500",  # Orange
        "#808080",  # Grey
        "#00FFFF",  # Cyan
        "#FF00FF",  # Magenta
        "#00FF00",  # Lime
        "#000080",  # Navy
        "#800000",  # Maroon
        "#808000",  # Olive
        "#008080",  # Teal
        "#C0C0C0",  # Silver
        "#FF6347",  # Tomato
        "#4682B4",  # SteelBlue
        "#D2691E",  # Chocolate
        "#FF4500",  # OrangeRed
        "#8A2BE2",  # BlueViolet
        "#A52A2A",  # Brown
        "#DEB887",  # BurlyWood
        "#5F9EA0",  # CadetBlue
        "#7FFF00",  # Chartreuse
        "#D2691E",  # Chocolate
        "#FF7F50",  # Coral
        "#6495ED",  # CornflowerBlue
        "#FFF8DC",  # Cornsilk
        "#DC143C",  # Crimson
        "#00FFFF",  # Cyan
        "#00008B",  # DarkBlue
        "#008B8B",  # DarkCyan
        "#B8860B",  # DarkGoldenRod
        "#A9A9A9",  # DarkGray
        "#006400",  # DarkGreen
        "#BDB76B",  # DarkKhaki
        "#8B008B",  # DarkMagenta
        "#556B2F",  # DarkOliveGreen
        "#FF8C00",  # DarkOrange
    ]  # 40 colors

    colors_full = [x.upper() for x in colors_full]

    # make a new color list to match the liquids order.
    # Use the same color for the same liquid.
    # If the liquid is used more than once,
    # use the same color for all instances.
    colors = []
    for i, liquid in enumerate(liquids):
        if liquid not in liquids[:i]:
            colors.append(colors_full[i])
        else:
            colors.append(colors[liquids.index(liquid)])

    for liquid, des, color, v, loc_list in zip(
        liquids, descriptions, colors, l_volumes, l_locations
    ):
        liq = protocol.define_liquid(
            name=str(liquid), description=str(des), display_color=color
        )
        for loc in loc_list:
            loc.load_liquid(liquid=liq, volume=v)

    # label the wells in the sample_plate with 50 µL volume
    for k, well in enumerate(sample_plate_wells):
        liq = protocol.define_liquid(
            f"RNA sample x{k+1}",
            f"RNA sample x{k+1}, 50 µL",
            display_color=colors_full[len(colors) + 1],
        )
        well.load_liquid(liquid=liq, volume=5)

    # label the reserved wells in the final_plate_p1 with 0 volume
    for k, well in enumerate(final_plate_wells):
        liq = protocol.define_liquid(
            f"Fragmented RNA x{k+1}",
            "Reserve for Fragmented RNA product, Part 1, no liquid when starting",
            display_color=colors_full[len(colors) + 2],
        )
        well.load_liquid(liquid=liq, volume=0)
    # ------------------------------------------------------------------ #
    #                       custom functions, basic                      #
    # ------------------------------------------------------------------ #

    # ------------------------- nonlocal params ------------------------ #
    pause_time = 1

    tip_count_200 = 0
    tiporder_200 = [t1 for i in range(len(tips200)) for t1 in tips200[i].rows()[0]]

    tip_count_50 = 0
    tiporder_50 = [t2 for i in range(len(tips50)) for t2 in tips50[i].rows()[0]]

    def pause_attention(msg="pause", flash=False):
        """pause the robot, flashes the light, and display a message"""
        nonlocal pause_time
        if flash:
            protocol.set_rail_lights(False)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(True)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(False)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(True)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(False)
            protocol.delay(seconds=0.25)
            protocol.set_rail_lights(True)
        protocol.comment(f"\n\n\nPAUSE X{pause_time}")
        protocol.pause(msg)
        pause_time += 1

    # -------------------------- tip handling -------------------------- #
    def pick_up_200():
        nonlocal tip_count_200, tiporder_200, tips200, tips200_extra
        if tip_count_200 == len(tiporder_200):
            if len(tips200_extra) == 0:
                pause_attention("Replace empty 200ul filter tips")
                tiporder_200.clear()
                for a in slots_200:
                    del protocol.deck[a]
                for c in expansion_slots_200:
                    del protocol.deck[c]
                tips200 = [
                    protocol.load_labware("opentrons_flex_96_filtertiprack_200ul", slot)
                    for slot in slots_200
                ]
                for a in range(len(tips200)):
                    tiporder_200 += tips200[a].rows()[0]
                tips200_extra = [
                    protocol.load_labware("opentrons_flex_96_filtertiprack_200ul", slot)
                    for slot in slots_200_extra
                ]
            else:  # still have extra 200 ul tips on the expansion slots
                tiporder_200.clear()
                if len(tips200_extra) == 4:
                    protocol.comment(f"\n---move out 1: from {slots_200[0]} to {chute}")
                    protocol.move_labware(tips200[0], chute, use_gripper=True)
                else:
                    protocol.comment(
                        f"\n---move out 1: from {slots_200[0]} to {slots_empty_expansion[0]}"
                    )
                    protocol.move_labware(
                        tips200[0], slots_empty_expansion[0], use_gripper=True
                    )
                for b in range(len(slots_200_extra)):
                    protocol.comment(
                        f"\n---move in {b+1}: from {slots_200_extra[b]} to {slots_200[b]}"
                    )
                    protocol.move_labware(
                        tips200_extra[b], slots_200[b], use_gripper=True
                    )
                    if b + 1 < len(slots_200_extra):
                        protocol.comment(
                            f"---move out {b+2}: from {slots_200[b+1]} to {slots_200_extra[b]}"
                        )
                        protocol.move_labware(
                            tips200[b + 1], slots_200_extra[b], use_gripper=True
                        )
                    tiporder_200 += tips200_extra[b].rows()[0]
                tips200_extra.clear()
            tip_count_200 = 0
        p1000m.pick_up_tip(tiporder_200[tip_count_200])
        tip_count_200 += 1

    def pick_up_50():
        nonlocal tip_count_50
        if tip_count_50 == len(tiporder_50):
            pause_attention(msg="\n\nReplace empty tipracks before resuming.")
            p50m.reset_tipracks()
            tip_count_50 = 0
        p50m.pick_up_tip(tiporder_50[tip_count_50])
        tip_count_50 += 1

    def tip_disposal(pip):
        if DRY_WATER_RUN:
            pip.return_tip()
        elif pip.has_tip:
            pip.drop_tip(chute)

    def slow_withdraw(pip, well, z, delay_seconds):
        pip.default_speed /= 40
        if delay_seconds > 0:
            protocol.delay(seconds=delay_seconds)
        pip.move_to(well.top(z))
        pip.default_speed *= 40

    # ------------------------------------------------------------------ #
    #                      custom function, advanced                     #
    # ------------------------------------------------------------------ #

    # -------------------------- custom mixing ------------------------- #
    def custom_mix(pip, mvol, mix_loc, mix_rep, blowout=False, low=False, high=False):
        # Adjustable low and high zoffset for bead mixing
        # The default is false for regular mixing
        if low:
            asp = mix_loc.bottom(low)
        else:
            asp = mix_loc.bottom(1.5)
        if high:
            disp = mix_loc.bottom(high)
        else:
            disp = mix_loc.bottom(2.5)

        # define mixing volume, and use the 2nd smallest value
        if pip == p50m:
            a = 2
        elif pip == p1000m:
            a = 5
        b = 0.8 * pip.tip_racks[0].wells()[0].max_volume
        c = mvol
        numbers = [a, b, c]
        vol = sorted(numbers)[1]

        for _ in range(mix_rep):
            pip.aspirate(vol, asp, rate=sample_rate)
            pip.dispense(vol, disp, push_out=0)
        if blowout:
            pip.flow_rate.blow_out /= 5
            slow_withdraw(pip, mix_loc, z=-3, delay_seconds=1)
            pip.blow_out(mix_loc.top(-3))
            pip.touch_tip(radius=0.75, v_offset=-3, speed=10)
            pip.flow_rate.blow_out *= 5

    # -------------------------- custom delay -------------------------- #
    def custom_delay(name, time):
        if DRY_WATER_RUN:
            for j in np.arange(time, 0, -time):
                msg = (
                    "There are "
                    + str(j)
                    + " seconds left in the "
                    + str(name)
                    + " step"
                )  # noqa
                protocol.delay(seconds=time, msg=msg)
        else:
            for j in np.arange(time, 0, -0.5):
                msg = (
                    "There are "
                    + str(j)
                    + " minutes left in the "
                    + str(name)
                    + " step"
                )  # noqa
                protocol.delay(minutes=0.5, msg=msg)

    # -------------------- well to list distribution ------------------- #
    def well_to_list(
        pip,
        well_list,
        dest_list,
        transfer_vol,
        liquid_rate,
        z_height,
        pre_mix=False,
        post_mix=False,
        post_mix_rep=0,
        post_mix_vol=0,
        tip_withdrawal=False,
        change_tip=False,
        index=0,
    ):
        for i, (well, d) in enumerate(zip(well_list, dest_list)):
            protocol.comment(
                f"\n---Column x{i+index+1}: {transfer_vol} µL from {well} to {d.parent}"
            )
            if change_tip:
                if pip == p1000m:
                    pick_up_200()
                elif pip == p50m:
                    pick_up_50()
            if pre_mix:
                custom_mix(pip, transfer_vol, well, 6)
            if pip == p50m and 1 <= transfer_vol <= 5:
                pip.configure_for_volume(transfer_vol)
            pip.aspirate(transfer_vol, well, rate=liquid_rate)
            pip.dispense(transfer_vol, d.bottom(z_height), rate=liquid_rate)
            if post_mix:
                # set a threshold for the post mix volume in case p50m in a low volume mode
                if pip == p50m and 1 <= transfer_vol <= 5:
                    post_mix_vol = 25 if post_mix_vol > 25 else post_mix_vol
                if post_mix_vol <= 10:
                    custom_mix(
                        pip, post_mix_vol, d, post_mix_rep, blowout=False, low=1, high=1
                    )
                else:
                    custom_mix(pip, post_mix_vol, d, post_mix_rep, blowout=False)
            if tip_withdrawal:
                slow_withdraw(pip, d, z=-3, delay_seconds=1)
            if pip == p50m and 1 <= transfer_vol <= 5:
                pip.configure_for_volume(50)
            if z_height < 0.5 * d.depth:
                pip.blow_out(d.top(-3))
            if change_tip:
                tip_disposal(pip)
            # else:
            #     if i == len(dest_list) - 1:
            #         tip_disposal(pip)

    # # ------------------------- Remove liquids ------------------------- #

    def remove(re_vol, re_s_list, mode="elution", change_tip=False):
        if mode not in [
            "beads_supernatant",
            "ethanol",
            "elution",
        ]:
            raise ValueError("A wrong mode is used")
        if mode == "elution":
            dest_list = final_plate_cols
        else:
            dest_list = waste_well * cols
        for i, (s, d) in enumerate(zip(re_s_list, dest_list)):
            protocol.comment(f"\n~~~Column X{i+1}: {mode.title()} Removal")
            protocol.comment(f"---Source: {str(s)}")
            protocol.comment(f"---Destination: {str(d)}")
            if mode == "elution":
                re_rate = elution_rate
                loc = d.bottom(1.5)
            elif mode == "ethanol":
                re_rate = ethanol_rate
                loc = d.top(-3)
            elif mode == "beads_supernatant":
                re_rate = sample_rate
                loc = d.top(-3)
            # define volume parameters
            if re_vol > 200:
                raise ValueError("A wrong volume is used")
            elif 180 < re_vol <= 200:
                pip = p1000m
                if not pip.has_tip:
                    pick_up_200()
                pip.aspirate(
                    50,
                    s.bottom().move(types.Point(x=0, y=0, z=0.5 * (s.depth))),
                    rate=0.1,
                )
                slow_withdraw(pip, s, z=3, delay_seconds=1)
                pip.dispense(50, loc, rate=re_rate, push_out=0)
                asp_vol = re_vol - 50
                extra_vol = 20
            elif 50 <= re_vol <= 180:
                pip = p1000m
                asp_vol = re_vol
                extra_vol = 20
            else:
                pip = p50m
                asp_vol = re_vol
                extra_vol = 5
            if asp_vol + extra_vol * 2 > pip.tip_racks[0].wells()[0].max_volume:
                extra_vol = (pip.tip_racks[0].wells()[0].max_volume - asp_vol) / 2
            disp_vol = asp_vol + extra_vol

            if not pip.has_tip:
                if pip == p1000m:
                    pick_up_200()
                else:
                    pick_up_50()
            # air gap for elution
            if mode == "elution":
                pip.move_to(s.top())
                pip.air_gap(extra_vol)
            # regular removal
            pip.aspirate(
                asp_vol,
                s.bottom().move(types.Point(x=0, y=0, z=0.7)),
                rate=0.1,
            )
            # extra removal for buffers and ethanol, air gap built up on the bottom
            if mode != "elution":
                pip.aspirate(
                    extra_vol,
                    s.bottom().move(types.Point(x=0, y=0, z=0.5)),
                    rate=0.05,
                )
            slow_withdraw(pip, s, z=-3, delay_seconds=1)
            pip.dispense(disp_vol, loc, rate=re_rate)  # no push out used
            protocol.delay(seconds=3)
            if mode == "elution":
                slow_withdraw(pip, d, z=-3, delay_seconds=0)
                pip.touch_tip(radius=0.75, v_offset=-3, speed=5)
            else:
                pip.blow_out(loc)
                pip.air_gap(extra_vol)
            if change_tip:
                tip_disposal(pip)
            else:
                if i == len(dest_list) - 1:
                    tip_disposal(pip)
            # if b2slots:
            #     protocol.move_labware(
            #         re_s_list[0].parent, transition_slot, use_gripper=True
            #     )
            protocol.comment("\n")

    # # ------------------------------ wash ------------------------------ #
    def wash(wash_vol, wash_s_list, addition_change_tip=False):
        height = sample_plate.wells()[0].depth
        if cols <= 8:
            s1 = [wash_s_list[0]] * cols
            dest_list_1 = sample_plate_cols
            s2 = None
            dest_list_2 = None
        else:
            s1 = [wash_s_list[0]] * 8
            dest_list_1 = sample_plate_cols[:8]
            s2 = [wash_s_list[1]] * (cols - 8)
            dest_list_2 = sample_plate_cols[8:]
        protocol.comment("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        protocol.comment(" Ethanol Wash starts")
        protocol.comment("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        if not addition_change_tip:
            pick_up_200()
        well_to_list(
            p1000m,
            s1,
            dest_list_1,
            wash_vol,
            ethanol_rate,
            height + 0.5,
            pre_mix=False,
            post_mix=False,
            post_mix_rep=0,
            post_mix_vol=0,
            tip_withdrawal=False,
            change_tip=addition_change_tip,
            index=0,
        )
        if s2 is not None and dest_list_2 is not None:
            well_to_list(
                p1000m,
                s2,
                dest_list_2,
                wash_vol,
                ethanol_rate,
                height + 0.5,
                pre_mix=False,
                post_mix=False,
                post_mix_rep=0,
                post_mix_vol=0,
                tip_withdrawal=False,
                change_tip=addition_change_tip,
                index=8,
            )
        if not addition_change_tip:
            tip_disposal(p1000m)
        custom_delay("Ethanol Wash Incubation", time=0.5)
        protocol.comment("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        protocol.comment(" Ethanol Removal starts")
        protocol.comment("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        remove(
            wash_vol,
            sample_plate_cols,
            mode="ethanol",
            change_tip=True,
        )

    def elution(volume_in, volume_out):
        well_to_list(
            p50m,
            [te] * cols,
            sample_plate_cols,
            volume_in,
            buffer_rate,
            0.8,
            pre_mix=False,
            post_mix=True,
            post_mix_rep=10,
            post_mix_vol=0.8 * volume_in,
            tip_withdrawal=True,
            change_tip=True,
        )
        custom_delay("Elution Incubation", time=beads_incubation_time)
        protocol.move_labware(sample_plate, mag, use_gripper=True)
        custom_delay("Elution Binding", time=beads_binding_time)
        thermocycler.open_lid()
        remove(
            volume_out,
            sample_plate_cols,
            mode="elution",
            change_tip=True,
        )
        thermocycler.close_lid()

    # ------------------------------------------------------------------ #
    #                        protocol starts here                        #
    # ------------------------------------------------------------------ #
    protocol.comment("\n---Protocol starts---\n")

    protocol.comment(
        f"\n---1A9.1 to 1A9.5: Resuspension Beads({res_beads_vol} uL) Addition, mixing, and supernatant removal---\n"
    )
    well_to_list(
        p1000m,
        [res_beads] * cols,
        sample_plate_cols,
        res_beads_vol,
        beads_rate,
        2,
        pre_mix=True,
        post_mix=True,
        post_mix_rep=10,
        post_mix_vol=0.8 * (res_beads_vol + sample_vol),
        tip_withdrawal=True,
        change_tip=True,
    )
    custom_delay("Resuspension Beads Incubation", time=beads_incubation_time_long)
    protocol.move_labware(sample_plate, mag, use_gripper=True)
    custom_delay("Resuspension Beads Binding", time=beads_binding_time)
    remove(
        res_beads_vol + sample_vol,
        sample_plate_cols,
        mode="beads_supernatant",
        change_tip=True,
    )
    protocol.move_labware(sample_plate, transition_slot, use_gripper=True)

    protocol.comment(
        f"\n---1A9.6 to 1A9.9: TE Buffer ({te_1_vol} µL)---\n"
        f"---NEBNext Bead Reconstitution Buffer ({recon_buffer_vol} µL)---\n"
        "---Addition, mixing, and supernatant removal---\n"
    )
    well_to_list(
        p1000m,
        [te] * cols,
        sample_plate_cols,
        te_1_vol,
        buffer_rate,
        1.5,
        pre_mix=False,
        post_mix=True,
        post_mix_rep=10,
        post_mix_vol=0.8 * te_1_vol,
        tip_withdrawal=True,
        change_tip=True,
    )
    well_to_list(
        p1000m,
        [recon_buffer] * cols,
        sample_plate_cols,
        recon_buffer_vol,
        buffer_rate,
        1.5,
        pre_mix=False,
        post_mix=True,
        post_mix_rep=10,
        post_mix_vol=0.8 * (recon_buffer_vol + te_1_vol),
        tip_withdrawal=True,
        change_tip=True,
    )
    custom_delay("TE+Recon Buffer Beads Incubation", time=beads_incubation_time_long)
    protocol.move_labware(sample_plate, mag, use_gripper=True)
    custom_delay("TE+Recon Buffer Beads Binding", time=beads_binding_time)
    remove(
        te_1_vol + recon_buffer_vol,
        sample_plate_cols,
        mode="beads_supernatant",
        change_tip=True,
    )

    protocol.comment(
        "\n---1A9.10: 1st round 80% Ethanol(200 µL) wash, and supernatant removal, stay on the magnetic block---\n"
    )
    wash(ethanol_vol, ethanol_1, addition_change_tip=False)

    protocol.comment(
        "\n---1A9.11 to 1A9.12: 2nd round 80% Ethanol (200 µL) wash, and supernatant removal, air dry---\n"
        "---move to the top of the temperature module---\n"
    )
    wash(ethanol_vol, ethanol_2, addition_change_tip=False)
    custom_delay("air dry", time=air_dry_time)
    protocol.move_labware(sample_plate, transition_slot, use_gripper=True)

    protocol.comment(
        f"\n---1A9.13 to 1A9.14: Elution(TE), {te_2_in_vol} µL in, {te_2_out_vol} µL out---\n"
    )
    elution(te_2_in_vol, te_2_out_vol)

    pause_attention("The protocol is complete. Please remove the final plate.")
    thermocycler.open_lid()
    protocol.move_labware(final_plate, transition_slot, use_gripper=True)
    thermocycler.close_lid()
    thermocycler.deactivate_lid()
    thermocycler.deactivate_block()
    temp_mod.deactivate()
