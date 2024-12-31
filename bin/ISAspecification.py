import math
import numpy as np

# 定数
R = 287.05  # 空気のガス定数 (J/kg⋅K)
gamma = 1.4  # 空気の比熱比
T0 = 288.15  # 海面上の気温 (K)
P0 = 101325  # 海面上の大気圧 (Pa)
mu_0 = 1.7894e-5  # 海面上の動粘性係数 (kg/m⋅s)
g = 9.80665  # 重力加速度 (m/s²)
L = 0.0065  # 温度減率 (K/m)

def get_temperature(altitude, dT=0):
    """
    高度に基づいて気圧と気温を計算する（国際標準大気に基づく）。
    altitude: 高度 (m)
    return: (pressure, temperature) 気圧 (Pa), 気温 (K)
    """
    if altitude <= 11000:  # 対流圏（高度 11,000 m まで）
        temperature_K = T0 - L * altitude
    else:  # 成層圏（11,000 m から 20,000 m の間）
        temperature_K = 216.65  # 一定温度
    return temperature_K + dT


def get_pressure(altitude, dT=0):
    """
    高度に基づいて気圧と気温を計算する（国際標準大気に基づく）。
    altitude: 高度 (m)
    return: (pressure, temperature) 気圧 (Pa), 気温 (K)
    """
    if altitude <= 11000:  # 対流圏（高度 11,000 m まで）
        pressure = P0 * (1 - (L * altitude) / T0) ** (g / (R * L))
    else:  # 成層圏（11,000 m から 20,000 m の間）
        pressure = P0 * 0.22336 * math.exp(-g * (altitude - 11000) / (R * 216.65))
    return pressure

def mach_to_velocity(mach, altitude=0, dT=0):
    """
    Convert Mach number to velocity based on the International Standard Atmosphere (ISA).
    mach: Mach number
    altitude: Altitude (m)
    """

    # Get the temperature based on altitude
    temperature_K = get_temperature(altitude, dT)
    
    # Calculate the velocity of sound from the temperature (a = sqrt(gamma * R * T))
    sound_velocity = math.sqrt(gamma * R * temperature_K)
    
    # Calculate the velocity from the Mach number (v = M * a)
    velocity = mach * sound_velocity
    return velocity

def velocity_to_mach(velocity, dT=0, altitude=0):
    """
    Convert Mach number to velocity based on the International Standard Atmosphere (ISA).
    velocity: air speed (m/s)
    altitude: Altitude (m)
    """

    # Get the temperature based on altitude
    temperature_K = get_temperature(altitude, dT)
    
    # Calculate the velocity of sound from the temperature (a = sqrt(gamma * R * T))
    sound_velocity = np.sqrt(gamma * R * temperature_K)
    
    # Calculate the velocity from the Mach number (v = M * a)
    mach = velocity / sound_velocity
    return mach

def feet_to_meters(feet):
    """
    Convert feet to meters.
    feet: Input length in feet (float or int)
    return: Length in meters (float)
    """
    meters = feet * 0.3048  # Conversion factor from feet to meters
    return meters

def get_density(altitude, dT):
    """
    dT: 気温 (°C) - 高度での気温
    altitude: 高度 (m)
    """
    
    # 高度に基づいて気圧と気温を取得
    temperature_K = get_temperature(altitude, dT)
    pressure = get_pressure(altitude, dT)
    
    # 空気密度を計算
    density = pressure / (R * temperature_K)
    
    return density

def velocity_to_reynolds(velocity, length, dT=0, altitude=0):
    """
    高度、温度、速度、代表長に基づいたレイノルズ数を計算する関数。
    dT: 気温 (°C) - 高度での気温
    velocity: 速度 (m/s)
    length: 代表長 (m)
    altitude: 高度 (m)
    """
    
    # 高度に基づいて気圧と気温を取得
    temperature_K = get_temperature(altitude, dT)
    pressure = get_pressure(altitude, dT)
    
    # 空気密度を計算
    density = pressure / (R * temperature_K)
    
    # 動粘性係数 (Sutherlandの式に基づく)
    C1 = 120  # Sutherland定数
    mu = mu_0 * (temperature_K / T0) ** (3/2) * (T0 + C1) / (temperature_K + C1)
    
    # レイノルズ数の計算
    reynolds = density * velocity * length / mu

    return reynolds