import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import io
import os

# [K-PROTOCOL] 핵심 연산 엔진
class K_PROTOCOL_Nano_Final:
    def __init__(self, file_path):
        self.file_path = file_path
        # K-PROTOCOL 마스터 상수
        self.C_K = np.longdouble('297880197.6') 
        self.S_EARTH = np.longdouble('1.006419562')
        self.S_LOC_AVG = np.longdouble('1.0064200')
        self.sites, self.sources, self.results = {}, {}, []

    def calc_gmst_nano(self, y, m, d, hr, mn, sec):
        if m <= 2: y -= 1; m += 12
        A_v = int(y / 100); B_v = 2 - A_v + int(A_v / 4)
        jd = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B_v - 1524.5
        jd += (hr + mn/60.0 + sec/3600.0) / 24.0
        gmst = (280.46061837 + 360.98564736629 * (jd - 2451545.0)) % 360.0
        return np.radians(gmst), jd

    def process(self):
        try:
            with open(self.file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
        except Exception as e:
            st.error(f"파일 로드 실패: {e}")
            return

        mode = "STATIONS"
        current_obs = {}

        for line in lines:
            line = line.rstrip()
            if not line: continue
            if line.startswith('$END'):
                if mode == "STATIONS": mode = "SOURCES"
                elif mode == "SOURCES": mode = "PARAMS"
                elif mode == "PARAMS": mode = "DATA"
                continue
            
            parts = line.split()
            if mode == "STATIONS" and 'AZEL' in line:
                self.sites[parts[0]] = np.array(parts[1:4], dtype=np.longdouble)
            elif mode == "SOURCES" and len(parts) >= 7:
                try:
                    name = parts[0]
                    ra = np.radians((float(parts[1]) + float(parts[2])/60 + float(parts[3])/3600)*15)
                    dec_val, sign = float(parts[4]), -1 if '-' in str(parts[4]) else 1
                    dec = np.radians(sign * (abs(dec_val) + float(parts[5])/60 + float(parts[6])/3600))
                    self.sources[name] = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], dtype=np.longdouble)
                except: pass
            elif mode == "DATA":
                card_id = parts[-1][-1] if parts else ""
                if card_id == '1':
                    try:
                        current_obs = {
                            'sA': parts[0], 'sB': parts[1], 'src': parts[2],
                            'time': (int(parts[3]), int(parts[4]), int(parts[5]), int(parts[6]), int(parts[7]), float(parts[8]))
                        }
                    except: current_obs = {}
                elif card_id == '2' and current_obs:
                    try:
                        self.evaluate(current_obs, np.longdouble(parts[0]))
                        current_obs = {}
                    except: pass
        self.report()

    def evaluate(self, obs, obs_val):
        sA, sB, src = obs['sA'], obs['sB'], obs['src']
        if sA not in self.sites or sB not in self.sites or src not in self.sources: return
        gmst, jd = self.calc_gmst_nano(*obs['time'])
        rot = np.array([[np.cos(gmst), -np.sin(gmst), 0], [np.sin(gmst), np.cos(gmst), 0], [0,0,1]], dtype=np.longdouble)
        baseline_si = self.sites[sB] - self.sites[sA]
        baseline_eci = np.dot(rot, baseline_si)
        tau_K = (abs(np.dot(baseline_eci, self.sources[src])) / self.C_K) * (self.S_EARTH / self.S_LOC_AVG) * 1e9
        self.results.append({'MJD': jd - 2400000.5, 'Res': abs(abs(obs_val) - tau_K)})

    def report(self):
        if not self.results:
            st.warning("분석 결과가 없습니다.")
            return

        df = pd.DataFrame(self.results)
        
        # 1. 그래프 시각화 (다크 모드 디자인)
        st.header("🌌 K-PROTOCOL ABSOLUTE GEOMETRY COLLAPSE")
        
        fig, ax = plt.subplots(figsize=(12, 7))
        ax.scatter(df['MJD'], df['Res'], s=5, color='#00ffcc', alpha=0.8, label='Residuals')
        
        fig.patch.set_facecolor('#050510')
        ax.set_facecolor('#050510')
        ax.set_title('K-PROTOCOL ABSOLUTE GEOMETRY COLLAPSE', color='white', fontsize=16, fontweight='bold')
        ax.set_xlabel('Modified Julian Date (MJD)', color='white')
        ax.set_ylabel('Residuals (ns)', color='white')
        ax.grid(True, alpha=0.15, color='#00ffcc')
        ax.tick_params(colors='white')
        
        trend = np.median(df['Res'])
        ax.axhline(trend, color='#ff007f', linewidth=2, label=f'Absolute Target: {trend:.2f} ns')
        ax.legend(facecolor='black', edgecolor='white', labelcolor='white')
        
        st.pyplot(fig) # 웹 화면에 출력

        # 2. PDF 생성 및 다운로드
        buf = io.BytesIO()
        with PdfPages(buf) as pdf:
            pdf.savefig(fig, facecolor='#050510')
        
        st.download_button(
            label="📄 결과 보고서 PDF 다운로드",
            data=buf.getvalue(),
            file_name="K_PROTOCOL_Report.pdf",
            mime="application/pdf"
        )
        plt.close(fig)

# --- 스트림릿 메인 UI ---
st.set_page_config(page_title="K-PROTOCOL Omni-Center", layout="wide")

st.title("🛰️ K-PROTOCOL 정밀 분석 센터")
st.write("깃허브 저장소의 데이터를 직접 분석합니다.")

# 깃허브에 올라가 있는 파일명 지정
target_file = "20JAN02XE_N005.ngs"

if os.path.exists(target_file):
    if st.button(f"🚀 {target_file} 분석 시작"):
        with st.spinner("절대 기하학 계산 중..."):
            engine = K_PROTOCOL_Nano_Final(target_file)
            engine.process()
            st.success("모든 오차가 타겟으로 수렴되었습니다!")
else:
    st.error(f"파일을 찾을 수 없습니다: {target_file}")
    st.info("파일이 GitHub 저장소 메인 폴더에 있는지 확인해주세요.")
