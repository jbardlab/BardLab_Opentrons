from opentrons import protocol_api

# from opentrons import types
from opentrons.protocol_api import SINGLE, ALL
import math
import numpy as np

metadata = {
    "protocolName": "NEBNext UltraExpress® RNA Library Prep Kit_NEB #E3330S/L_Part3_V1",
    "author": "Opentrons",
    "description": "Part 3: End Prep of cDNA Library, Adaptor Ligation, and PCR Enrichment of Adaptor Ligated DNA",
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
    # parameters.add_float(
    #     variable_name="beads_rate",
    #     display_name="Beads Rate",
    #     description="Beads mixing rate",
    #     default=0.2,
    #     minimum=0,
    #     maximum=1.5,
    # )
    parameters.add_float(
        variable_name="buffer_rate",
        display_name="Buffer Rate",
        description="Buffer mixing rate",
        default=0.8,
        minimum=0,
        maximum=1.5,
    )
    # parameters.add_float(
    #     variable_name="ethanol_rate",
    #     display_name="Ethanol Rate",
    #     description="Ethanol mixing rate",
    #     default=0.8,
    #     minimum=0,
    #     maximum=1.5,
    # )
    parameters.add_float(
        variable_name="enzyme_rate",
        display_name="Enzyme Rate",
        description="Enzyme mixing rate",
        default=0.5,
        minimum=0,
        maximum=1.5,
    )
    parameters.add_int(
        variable_name="enrichment_cycles",
        display_name="Enrichment Cycles(10-14)",
        description="The number of cycles for the enrichment step",
        default=12,
        minimum=10,
        maximum=14,
        unit="x Cycles",
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
    # beads_rate = protocol.params.beads_rate
    buffer_rate = protocol.params.buffer_rate
    # ethanol_rate = protocol.params.ethanol_rate
    enzyme_rate = protocol.params.enzyme_rate
    enrichment_cycles = protocol.params.enrichment_cycles
    # ------------------------ volume parameters ----------------------- #
    sample_vol = 20
    end_prep_mm_vol = 4
    diluted_adapter_vol = 2
    ligation_mm_vol = 12
    user_ez_vol = 2
    mstc_primer_mm_vol = 60

    # -------------------- magnetic block parameters ------------------- #
    # beads_incubation_time = 2
    # beads_binding_time = 5
    # air_dry_time = 5

    # ------------------------------------------------------------------ #
    #                           Error Handling                           #
    # ------------------------------------------------------------------ #

    # 1. Check if the number of samples is within the range
    if num_samples > 96 or num_samples < 1:
        raise ValueError("The number of samples should be between 1 and 96")
    # 2. Any of the volumes should be equal to or less than 200 µL
    if (
        sample_vol > 200
        or end_prep_mm_vol > 200
        or diluted_adapter_vol > 200
        or ligation_mm_vol > 200
        or user_ez_vol > 200
        or mstc_primer_mm_vol > 200
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
    # Rule 1: Use 50 µL tips for the 2x TC prep and elution steps.
    # Rule 2: Other steps use 200 µL tips for better beads mixing and removal.

    # --------------------- step by step analysis --------------------- #

    # ------------------- calculate 200 ul tips first ------------------ #
    # 5. steps 1A8.1 to 1A8.3: 50 µL NEBnext MSTC High Yield MM + 10 µL Index Primer combo (60 µL) in,
    # cols x columns: pre and post mixing, no liquid removal
    # TC steps: Lid to 105°C, 1 cycle of 98°C for 30 sec, 12 cycles of 98°C for 10 sec, 65°C for 75 sec, 1 cycle of 65°C for 5 min, hold at 4°C
    required_cols_200_1 = cols

    # ------------------- calculate 50 ul tips first ------------------- #

    # 1. steps 1A6.1 to 1A6.4: 2.5 µL End Prep buffer + 1.5 µL End Prep Enzyme mix (4 µL) in,
    # cols x columns: pre and post mixing, no liquid removal
    # TC steps: Lid to 75°C, 5min at 20°C, 10 min at 65°C, hold at 4°C
    required_cols_50_1 = cols

    # 2. steps 1A7.1 to 1A7.4: 2 µL diluted Adapter (2 µL) + 12 µL Ligation MM (separately) in,
    # cols x columns: pre and post mixing, no liquid removal
    # TC steps: Lid to 40°C, 15 min at 20°C, hold at 4°C
    required_cols_50_2 = cols + cols

    # 3. steps 1A7.5 to 1A7.6: 2 µL USER Enzyme (2 µL) in
    # cols x columns: pre and post mixing, no liquid removal
    # TC steps: Lid to 45°C, 5 min at 37°C, hold at 4°C
    required_cols_50_3 = cols

    # ---------------------- calculate 200 ul tips --------------------- #
    required_cols_200 = required_cols_200_1
    required_slots_num_200 = math.ceil(required_cols_200 / 12)
    # Assign 200 µl tips to C2
    if required_slots_num_200 > 1:
        raise ValueError(
            "200 µl tips exceed available space in C2"
        )  # D2 can only hold one rack
    slots_200 = ["D2"]

    # ---------------------- calculate 50 ul tips ---------------------- #
    required_cols_50 = required_cols_50_1 + required_cols_50_2 + required_cols_50_3
    required_slots_num_50 = math.ceil(required_cols_50 / 12)
    # ------------------- assign slots for 50ul tips ------------------ #
    regular_slots_50 = ["A2", "A3"]  # Regular slots available for 50 µl tips
    expansion_slots_50 = ["A4", "B4", "C4", "D4"]  # Expansion slots for extra racks

    difference = required_slots_num_50 - len(regular_slots_50)

    if difference <= 0:
        # All required racks fit in the regular slots
        slots_50 = regular_slots_50[:required_slots_num_50]
        slots_50_extra = []  # No overflow
        slots_empty_expansion = (
            expansion_slots_50  # All expansion slots remain available
        )
    else:
        # Use all available regular slots, overflow into expansion slots
        slots_50 = regular_slots_50  # Fill all regular slots
        slots_50_extra = expansion_slots_50[:difference]  # Use required expansion slots
        slots_empty_expansion = expansion_slots_50[
            difference:
        ]  # Leftover expansion slots

    print(f"Number of samples: {num_samples} x samples")
    print(f"50ul tips number: {required_slots_num_50} x slots")
    print(f"50ul tips columns: {int(cols*2)}")
    print(f"50ul tips slots: {slots_50}")
    print(f"50ul tips extra slots: {slots_50_extra}")
    print(f"200ul tips number: {required_slots_num_200} x slots")
    print(f"200ul tips columns: {int(required_cols_200)}")
    print(f"200ul tips slots: {slots_200}")
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
    # ------------------------------------------------------------------ #
    #                               Modules                              #
    # ------------------------------------------------------------------ #
    thermocycler = protocol.load_module("thermocycler module gen2")
    thermocycler.set_block_temperature(temperature=4)
    thermocycler.open_lid()
    mag = protocol.load_module("magneticBlockV1", "D1")
    temp_mod = protocol.load_module("temperature module gen2", "C1")
    temp_adapter = temp_mod.load_adapter(
        "opentrons_96_well_aluminum_block"
    )  # use the 96 PCR plate for better support low-volume reagents
    temp_mod.set_temperature(celsius=4)
    chute = protocol.load_waste_chute()
    # ------------------------------------------------------------------ #
    #                               labware                              #
    # ------------------------------------------------------------------ #
    reagent_plate = temp_adapter.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt", label="Reagent Plate"
    )
    sample_plate = thermocycler.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt", "Sample Plate"
    )
    index_plate = protocol.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt",
        "C2",
        "index plate",
    )
    # waste_res = protocol.load_labware("nest_1_reservoir_195ml", "D2", "Waste Reservoir")
    transition_slot = mag
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
    end_prep_mm = reagent_plate.rows()[0][0]  # A1
    end_prep_mm_list = reagent_plate.columns()[0]
    diluted_adapter = reagent_plate.rows()[0][1]  # A2
    diluted_adapter_list = reagent_plate.columns()[1]
    ligation_mm = reagent_plate.rows()[0][2]  # A3
    ligation_mm_list = reagent_plate.columns()[2]
    user_ez = reagent_plate.rows()[0][3]  # A4
    user_ez_list = reagent_plate.columns()[3]
    # mstc primer mm is assigned according to the number of columns
    #
    # if cols <= 3:
    #     mstc_primer_mm = [reagent_plate.rows()[0][4]]  # A5
    #     mstc_primer_mm_list = reagent_plate.columns()[4]
    # elif 4 <= cols <= 6:
    #     mstc_primer_mm = reagent_plate.rows()[0][4:6]  # A5, A6
    #     mstc_primer_mm_list = reagent_plate.columns()[4] + reagent_plate.columns()[5]
    # elif 7 <= cols <= 9:
    #     mstc_primer_mm = reagent_plate.rows()[0][4:7]  # A5, A6, A7
    #     mstc_primer_mm_list = (
    #         reagent_plate.columns()[4]
    #         + reagent_plate.columns()[5]
    #         + reagent_plate.columns()[6]
    #     )
    # elif 10 <= cols <= 12:
    #     mstc_primer_mm = reagent_plate.rows()[0][4:8]  # A5, A6, A7, A8
    #     mstc_primer_mm_list = (
    #         reagent_plate.columns()[4]
    #         + reagent_plate.columns()[5]
    #         + reagent_plate.columns()[6]
    #         + reagent_plate.columns()[7]
    # )
    mstc_primer_mm = index_plate.rows()[0][:cols]
    mstc_primer_mm_list = index_plate.wells()[:num_samples]
    # ------------------------- waste reservoir ------------------------ #
    # waste_well = [waste_res.wells()[0]]

    # summarize reagent plate map in a list
    reagent_plate_list = [
        "Column 1: End Prep Master Mix (2.5 µL End Prep buffer + 1.5 µL End Prep Enzyme mix per sample)",
        "Column 2: Diluted Adapter (2 µL per sample)",
        "Column 3: Ligation Master Mix (12 µL per sample)",
        "Column 4: USER Enzyme (2 µL per sample)",
    ]
    print(reagent_plate_list)
    # # # ------------------------------------------------------------------ #
    # # #                               liquids                              #
    # # # ------------------------------------------------------------------ #
    l_locations = [
        end_prep_mm_list,
        diluted_adapter_list,
        ligation_mm_list,
        user_ez_list,
        mstc_primer_mm_list,
    ]
    l_volumes = [
        end_prep_mm_vol * cols * 1.2,
        diluted_adapter_vol * cols * 1.2,
        ligation_mm_vol * cols * 1.2,
        user_ez_vol * cols * 1.2,
        mstc_primer_mm_vol * 1.1,
    ]
    liquids = [
        "End Prep Master Mix",
        "Diluted Adapter",
        "Ligation Master Mix",
        "USER Enzyme",
        "NEBNext MSTC High Yield MM + Index Primer Mix",
    ]
    descriptions = [
        "NEBNext Ultra II End Prep Master Mix, 4 µL per sample",
        "Diluted Adapter, 2 µL per sample",
        "NEBNext Ultra II Ligation Master Mix, 12 µL per sample",
        "USER Enzyme, 2 µL per sample",
        "NEBNext MSTC High Yield MM + Index Primer Mix,\n"
        "each well contains 10% extra liquid,66 µL = 55 µL NEBnext MSTC High Yield MM + 11 µL Index Primer\n"
        "When in use, 60 µL per sample (different index for each individual sample)\n",
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
            f"cDNA sample x{k+1} from Part 2, 20 µL",
            f"cDNA sample x{k+1} from Part 2, 20 µL",
            display_color=colors_full[len(colors) + 1],
        )
        well.load_liquid(liquid=liq, volume=20)

    # # # ------------------------------------------------------------------ #
    # # #                       custom functions, basic                      #
    # # # ------------------------------------------------------------------ #

    # # # ------------------------- nonlocal params ------------------------ #
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

    # # -------------------------- tip handling -------------------------- #
    def pick_up_50():
        nonlocal tip_count_50, tiporder_50, tips50, tips50_extra
        if tip_count_50 == len(tiporder_50):
            if len(tips50_extra) == 0:
                pause_attention("Replace empty 50ul filter tips")
                tiporder_50.clear()
                for a in slots_50:
                    del protocol.deck[a]
                for c in expansion_slots_50:
                    del protocol.deck[c]
                tips50 = [
                    protocol.load_labware("opentrons_flex_96_filtertiprack_50ul", slot)
                    for slot in slots_50
                ]
                for a in range(len(tips50)):
                    tiporder_50 += tips50[a].rows()[0]
                tips50_extra = [
                    protocol.load_labware("opentrons_flex_96_filtertiprack_50ul", slot)
                    for slot in slots_50_extra
                ]
            else:  # still have extra 200 ul tips on the expansion slots
                tiporder_50.clear()
                if len(tips50_extra) == 3:
                    protocol.comment(f"\n---move out 1: from {slots_50[0]} to {chute}")
                    protocol.move_labware(tips50[0], chute, use_gripper=True)
                else:
                    protocol.comment(
                        f"\n---move out 1: from {slots_50[0]} to {slots_empty_expansion[0]}"
                    )
                    protocol.move_labware(
                        tips50[0], slots_empty_expansion[0], use_gripper=True
                    )
                for b in range(len(slots_50_extra)):
                    protocol.comment(
                        f"\n---move in {b+1}: from {slots_50_extra[b]} to {slots_50[b]}"
                    )
                    protocol.move_labware(
                        tips50_extra[b], slots_50[b], use_gripper=True
                    )
                    if b + 1 < len(slots_50_extra):
                        protocol.comment(
                            f"---move out {b+2}: from {slots_50[b+1]} to {slots_50_extra[b]}"
                        )
                        protocol.move_labware(
                            tips50[b + 1], slots_50_extra[b], use_gripper=True
                        )
                    tiporder_50 += tips50_extra[b].rows()[0]
                tips50_extra.clear()
            tip_count_50 = 0
        p50m.pick_up_tip(tiporder_50[tip_count_50])
        tip_count_50 += 1

    def pick_up_200():
        nonlocal tip_count_200
        if tip_count_200 == len(tiporder_200):
            pause_attention(msg="\n\nReplace empty tipracks before resuming.")
            p1000m.reset_tipracks()
            tip_count_200 = 0
        p1000m.pick_up_tip(tiporder_200[tip_count_200])
        tip_count_200 += 1

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
                slow_withdraw(pip, d, z=-3, delay_seconds=0)
            if pip == p50m and 1 <= transfer_vol <= 5:
                pip.configure_for_volume(50)
            if z_height < 0.5 * d.depth:
                pip.blow_out(d.top(-3))
            if change_tip:
                tip_disposal(pip)
            # else:
            #     if i == len(dest_list) - 1:
            #         tip_disposal(pip)

    # # ---------------------- master mix preparation --------------------- #
    def mm_prep(mm_vol, post_mix_vol, well_list, column_list, change_partial_tip=False):
        nonlocal tip_count_50, tip_count_200
        # 1. column-wise distribution of MM
        pipette = p50m if mm_vol <= 50 else p1000m
        well_to_list(
            pipette,
            well_list,
            sample_plate_cols_m,
            mm_vol,
            enzyme_rate,
            1.5,
            pre_mix=True,
            post_mix=True,
            post_mix_rep=20,
            post_mix_vol=post_mix_vol,
            tip_withdrawal=False,
            change_tip=True,
        )
        # if number of samples is not a multiple of 8, the last column will be handled by the single tip
        if cols_s:
            if pipette == p50m:
                if tip_count_50 > 12:
                    column_index = tip_count_50 % 12
                else:
                    column_index = tip_count_50
                tips50_available = tiporder_50[tip_count_50].parent.columns()[
                    column_index
                ][::-1]
                if 1 <= mm_vol <= 5:
                    p50m.configure_for_volume(mm_vol)
                p50m.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=tips50)
                for i, (s, d) in enumerate(zip(column_list, sample_plate_cols_s)):
                    protocol.comment(
                        f"\n---Partial Column x{cols_m+1}, well x{i+1}: {mm_vol} µL"
                    )
                    if not p50m.has_tip:
                        p50m.pick_up_tip(location=tips50_available[i])
                    custom_mix(p50m, mm_vol, s, 10)
                    p50m.aspirate(mm_vol, s, rate=buffer_rate)
                    p50m.dispense(mm_vol, d.bottom(1.5), rate=buffer_rate, push_out=0)
                    if 1 <= mm_vol <= 5:
                        post_mix_vol = 25 if post_mix_vol > 25 else post_mix_vol
                    custom_mix(p50m, post_mix_vol, d, 10)
                    slow_withdraw(p50m, d, z=-3, delay_seconds=0)
                    p50m.blow_out(d.top(-3))
                    if change_partial_tip:
                        if not DRY_WATER_RUN:
                            tip_disposal(p50m)
                        else:
                            p50m.drop_tip(chute)
                    else:
                        if i == len(sample_plate_cols_s) - 1:
                            if not DRY_WATER_RUN:
                                tip_disposal(p50m)
                            else:
                                p50m.drop_tip(chute)
                if 1 <= mm_vol <= 5:
                    p50m.configure_for_volume(50)
                p50m.configure_nozzle_layout(style=ALL, tip_racks=tips50)
                tip_count_50 += 1
            elif pipette == p1000m:
                column_index = tip_count_200
                tips200_available = tiporder_200[tip_count_200].parent.columns()[
                    column_index
                ][::-1]
                for i, (s, d) in enumerate(zip(column_list, sample_plate_cols_s)):
                    protocol.comment(
                        f"\n---Partial Column x{cols_m+1}, well x{i+1}: {mm_vol} µL"
                    )
                    if not p1000m.has_tip:
                        p1000m.pick_up_tip(location=tips200_available[i])
                    custom_mix(p1000m, mm_vol, s, 10)
                    p1000m.aspirate(mm_vol, s, rate=buffer_rate)
                    p1000m.dispense(mm_vol, d.bottom(1.5), rate=buffer_rate, push_out=0)
                    custom_mix(p1000m, post_mix_vol, d, 10)
                    slow_withdraw(p1000m, d, z=-3, delay_seconds=0)
                    p1000m.blow_out(d.top(-3))
                    if change_partial_tip:
                        if not DRY_WATER_RUN:
                            tip_disposal(p1000m)
                        else:
                            p1000m.drop_tip(chute)
                    else:
                        if i == len(sample_plate_cols_s) - 1:
                            if not DRY_WATER_RUN:
                                tip_disposal(p1000m)
                            else:
                                p1000m.drop_tip(chute)
                p1000m.configure_nozzle_layout(style=ALL, tip_racks=tips200)
                tip_count_200 += 1

    # -------------------------- thermocycler steps --------------------- #
    def tc_steps(
        tc_vol,
        lid_temp,
        t1,
        time1,
        t2,
        time2,
        t3,
        time3,
        h_temp,
        custom_message=False,
    ):
        nonlocal tip_count_50, pause_time
        protocol.comment("\n---Thermocycler incubation---")
        if not DRY_WATER_RUN:
            thermocycler.set_lid_temperature(lid_temp)
        thermocycler.close_lid()
        if not DRY_WATER_RUN:
            thermocycler.set_block_temperature(
                t1, hold_time_minutes=time1, block_max_volume=tc_vol
            )
            if t2 and time2:
                thermocycler.set_block_temperature(
                    t2, hold_time_minutes=time2, block_max_volume=tc_vol
                )
            if t3 and time3:
                thermocycler.set_block_temperature(
                    t3, hold_time_minutes=time3, block_max_volume=tc_vol
                )
            thermocycler.set_block_temperature(h_temp)
        else:
            custom_delay("Dry run delay", time=1)
        thermocycler.open_lid()
        thermocycler.deactivate_lid()

    def tc_profile(lid_temp, profile, reps, reaction_vol, h_temp):
        protocol.comment(f"\n---Thermocycler Profile: {profile}---")
        thermocycler.set_lid_temperature(lid_temp)
        thermocycler.close_lid()
        thermocycler.set_block_temperature(
            98, hold_time_minutes=0.5, block_max_volume=reaction_vol
        )
        thermocycler.execute_profile(
            steps=profile, repetitions=reps, block_max_volume=reaction_vol
        )
        thermocycler.set_block_temperature(
            65, hold_time_minutes=5, block_max_volume=reaction_vol
        )
        thermocycler.set_block_temperature(h_temp)

    # # ------------------------------------------------------------------ #
    # #                        protocol starts here                        #
    # # ------------------------------------------------------------------ #
    protocol.comment("\n---Protocol starts---\n")

    protocol.comment("\n---1A6.1 to 1A6.4: End Prep(4 µL, low volume) and TC---\n")
    mm_prep(
        end_prep_mm_vol,
        (sample_vol + end_prep_mm_vol) * 0.8,
        [end_prep_mm] * cols,
        end_prep_mm_list,
        change_partial_tip=True,
    )
    tc_steps(
        tc_vol=end_prep_mm_vol + sample_vol,
        lid_temp=75,
        t1=20,
        time1=5,
        t2=65,
        time2=10,
        t3=None,
        time3=None,
        h_temp=4,
    )

    protocol.comment(
        "\n---1A7.1 to 1A7.4: Diluted Adapter(2 µL, low volume), Ligation MM(12 µL) and TC---\n"
    )
    mm_prep(
        diluted_adapter_vol,
        (diluted_adapter_vol + sample_vol + end_prep_mm_vol) * 0.8,
        [diluted_adapter] * cols,
        diluted_adapter_list,
        change_partial_tip=True,
    )
    mm_prep(
        ligation_mm_vol,
        (ligation_mm_vol + diluted_adapter_vol + sample_vol + end_prep_mm_vol) * 0.8,
        [ligation_mm] * cols,
        ligation_mm_list,
        change_partial_tip=True,
    )
    tc_steps(
        tc_vol=diluted_adapter_vol + ligation_mm_vol + sample_vol + end_prep_mm_vol,
        lid_temp=40,
        t1=20,
        time1=15,
        t2=None,
        time2=None,
        t3=None,
        time3=None,
        h_temp=4,
    )

    protocol.comment("\n---1A7.5 to 1A7.6: USER Enzyme(2 µL, low volume) and TC---\n")
    mm_prep(
        user_ez_vol,
        (
            user_ez_vol
            + ligation_mm_vol
            + diluted_adapter_vol
            + sample_vol
            + end_prep_mm_vol
        )
        * 0.8,
        [user_ez] * cols,
        user_ez_list,
        change_partial_tip=True,
    )
    tc_steps(
        tc_vol=user_ez_vol
        + ligation_mm_vol
        + diluted_adapter_vol
        + sample_vol
        + end_prep_mm_vol,
        lid_temp=45,
        t1=37,
        time1=5,
        t2=None,
        time2=None,
        t3=None,
        time3=None,
        h_temp=4,
    )

    protocol.comment("\n---1A7.7 to 1A7.9: NEBNext MSTC Primer Mix(60 µL) and TC---\n")
    pause_attention(
        "Freshly make 66 µL = 55 µL NEBnext MSTC High Yield MM + 11 µL Index Primer in each well of the index plate\n"
        "After resuming the protocol, The protocol will allow you move away the regular reagent plate first"
        "Then, it will guide you to transfer the updated index plate on the top of aluminum block from the temperature module"
    )
    protocol.move_labware(reagent_plate, protocol_api.OFF_DECK, use_gripper=False)
    protocol.move_labware(index_plate, temp_adapter, use_gripper=False)
    # mm_prep(
    #     mstc_primer_mm_vol,
    #     (
    #         mstc_primer_mm_vol
    #         + user_ez_vol
    #         + ligation_mm_vol
    #         + diluted_adapter_vol
    #         + sample_vol
    #         + end_prep_mm_vol
    #     )
    #     * 0.8,
    #     mstc_primer_mm * cols,
    #     mstc_primer_mm_list,
    #     change_partial_tip=True,
    # )
    well_to_list(
        p1000m,
        mstc_primer_mm,
        sample_plate_cols,
        mstc_primer_mm_vol,
        enzyme_rate,
        1.5,
        pre_mix=True,
        post_mix=True,
        post_mix_rep=20,
        post_mix_vol=(
            mstc_primer_mm_vol
            + user_ez_vol
            + ligation_mm_vol
            + diluted_adapter_vol
            + sample_vol
            + end_prep_mm_vol
        )
        * 0.8,
        tip_withdrawal=False,
        change_tip=True,
    )
    if not DRY_WATER_RUN:
        enrichment_profile = [
            {"temperature": 98, "hold_time_seconds": 10},
            {"temperature": 65, "hold_time_seconds": 75},
        ]
        tc_profile(
            lid_temp=105,
            profile=enrichment_profile,
            reps=enrichment_cycles,
            reaction_vol=mstc_primer_mm_vol
            + user_ez_vol
            + ligation_mm_vol
            + diluted_adapter_vol
            + sample_vol
            + end_prep_mm_vol,
            h_temp=4,
        )
    else:
        thermocycler.close_lid()
    pause_attention(
        "The protocol has been completed. Please remove the reagent plate on the Temperature module.\n"
        "Temperature module will be deactivated shortly to reach out to the room temperature.\n"
        "The robot will move the sample plate from the thermocycler to slot Magnetic Block\n"
        "The user need to put the sample plate on the temperature module after the last step.\n"
        "Place a new 96-well plate on the thermocycler module for the next step.\n"
        "Prepare all Room Temperature reagents for the next step.",
    )
    thermocycler.open_lid()
    protocol.move_labware(sample_plate, transition_slot, use_gripper=True)
    temp_mod.deactivate()
